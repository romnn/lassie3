"""Command line entry point for the baseline comparison."""

from __future__ import annotations

import argparse
import datetime as dt
import json
import logging
import sys
from dataclasses import replace
from pathlib import Path

from lassie3.settings import DEFAULT_DAY, RunSettings


def _settings(args: argparse.Namespace) -> RunSettings:
    settings = RunSettings(day=args.day, hours=args.hours, tag=args.tag)
    if args.lassie_spacing:
        settings = replace(settings, grid=replace(settings.grid, lassie_spacing=args.lassie_spacing))
    if args.velocity_model:
        settings = replace(settings, velocity_model=args.velocity_model)
    if args.blacklist:
        settings = replace(settings, stations_blacklist=tuple(args.blacklist))
    if getattr(args, "lassie_threshold", None) is not None:
        settings = replace(settings, lassie_detector_threshold=args.lassie_threshold)
    if getattr(args, "figures", False):
        settings = replace(settings, lassie_save_figures=True)
    if getattr(args, "lassie_mad_multiple", None) is not None:
        settings = replace(settings, lassie_mad_multiple=args.lassie_mad_multiple)
    if getattr(args, "qseek_threshold", None) is not None:
        settings = replace(settings, qseek_detection_threshold=args.qseek_threshold)
    if args.station_terms:
        settings = replace(settings, station_terms=args.station_terms)
    return settings


def _day(value: str) -> dt.date:
    return dt.date.fromisoformat(value)


def _threshold(keyword: str):
    """Parse a threshold that may be an absolute level or an adaptive keyword."""

    def parse(value: str) -> str | float:
        return value if value == keyword else float(value)

    return parse


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="lassie3",
        description="Baseline comparison of Lassie v1 and Qseek on ICDP Eger data.",
    )
    parser.add_argument("--day", type=_day, default=DEFAULT_DAY,
                        help=f"UTC day to process (default: {DEFAULT_DAY})")
    parser.add_argument("--hours", type=float, default=24.0,
                        help="length of the processed window from 00:00 UTC (default: 24)")
    parser.add_argument("--velocity-model", type=Path,
                        help="1D model in Pyrocko nd format, shared by both detectors")
    parser.add_argument("--blacklist", nargs="*", default=None, metavar="NET.STA",
                        help="stations to exclude from both runs")
    parser.add_argument("--tag", default="",
                        help="suffix for a variant run so it does not overwrite the baseline")
    parser.add_argument("--lassie-spacing", type=float, metavar="METRES",
                        help="Lassie grid node spacing (default from settings, 2000)")
    parser.add_argument("--station-terms", type=Path, metavar="TERMS_JSON",
                        help="terms.json from `ssst-terms`; both detectors then add the "
                             "source-specific station terms to their travel times")
    parser.add_argument("-v", "--verbose", action="store_true")

    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("prepare", help="stage waveforms, stations and velocity model")
    sub.add_parser("doctor", help="repair the Qseek/PyTorch OpenMP conflict")
    sub.add_parser("fetch-catalog", help="download the EMSC/BGR reference catalog")

    run_lassie = sub.add_parser("run-lassie", help="run the Lassie v1 search")
    run_lassie.add_argument("--lassie-threshold", type=_threshold("auto"),
                            help='absolute level on the stacked image function, or '
                                 '"auto" to place it above the run\'s own noise floor')
    run_lassie.add_argument("--lassie-mad-multiple", type=float,
                            help="MADs above the median at which \"auto\" sets the threshold")
    run_lassie.add_argument("--figures", action="store_true",
                            help="save a figure per detection (slow above a few hundred)")
    run_lassie.add_argument("--force", action="store_true", help="overwrite an existing run")

    run_qseek = sub.add_parser("run-qseek", help="run the Qseek search")
    run_qseek.add_argument("--qseek-threshold", type=_threshold("MAD"),
                           help='semblance threshold, or "MAD" for the adaptive default')
    run_qseek.add_argument("--force", action="store_true", help="overwrite an existing run")

    plot_cmd = sub.add_parser("plot", help="render the comparison figure")
    plot_cmd.add_argument("--lassie-min-mad", type=float,
                          help="show only Lassie detections above median + N MAD of the "
                               "run's own image function")
    plot_cmd.add_argument("--qseek-min-semblance", type=float,
                          help="show only Qseek detections at or above this semblance")
    plot_cmd.add_argument("--suffix", default="",
                          help='appended to the figure name, e.g. "-strict"')
    plot_cmd.add_argument("--no-jitter", action="store_true",
                          help="draw detections exactly on their grid nodes instead of "
                               "jittered within the node cell")
    sub.add_parser("summary", help="print catalog statistics as JSON")
    sub.add_parser("review", help="adversarial checks of the results against the reference catalog")
    sweep_cmd = sub.add_parser("sweep", help="Lassie runs at several MAD levels, tabulated with Qseek post-filters")
    sweep_cmd.add_argument("--mad", type=float, nargs="+", default=[20.0, 40.0, 60.0],
                           help="MAD multiples to run (default: 20 40 60)")
    sweep_cmd.add_argument("--no-run", action="store_true",
                           help="only tabulate tagged runs that already exist")

    terms_cmd = sub.add_parser("ssst-terms", help="derive source-specific station terms from the "
                                                  "reference events of the other days")
    terms_cmd.add_argument("--radius", type=float, default=3000.0, metavar="METRES",
                           help="hypocentral separation beyond which a reference event stops "
                                "informing a node (default 3000)")
    terms_cmd.add_argument("--prior-weight", type=float, default=5.0,
                           help="weight of the static term, in reference-event equivalents (default 5)")
    terms_cmd.add_argument("--outlier-level", type=float, default=6.0,
                           help="reject residuals beyond this many robust standard deviations (default 6)")
    terms_cmd.add_argument("--min-probability", type=float, default=0.3,
                           help="minimum PhaseNet probability for a pick (default 0.3)")
    terms_cmd.add_argument("--train-days", type=_day, nargs="*", default=None, metavar="YYYY-MM-DD",
                           help="days to learn from (default: every day with data except DAY)")
    terms_cmd.add_argument("--limit", type=int,
                           help="use only the first N reference events, for a quick check")
    terms_cmd.add_argument("--force", action="store_true", help="re-derive even if terms.json exists")
    sub.add_parser("compare-ssst", help="four-way table and figure: each detector with and "
                                        "without station terms")

    args = parser.parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)-7s %(name)s: %(message)s",
    )

    settings = _settings(args)

    match args.command:
        case "prepare":
            from lassie3.prepare import prepare

            prepare(settings)
        case "doctor":
            from lassie3.doctor import doctor

            doctor()
        case "fetch-catalog":
            from lassie3.reference import fetch, load

            fetch(settings, force=True)
            print(f"{len(load(settings))} reference events for {settings.day_str}")
        case "run-lassie":
            from lassie3.run_lassie import run

            print(f"Lassie run directory: {run(settings, force=args.force)}")
        case "run-qseek":
            from lassie3.doctor import unify_openmp
            from lassie3.run_qseek import run

            unify_openmp()

            print(f"Qseek run directory: {run(settings, force=args.force)}")
        case "plot":
            from lassie3.plot import plot

            lassie_min = None
            if args.lassie_min_mad is not None:
                from lassie3.run_lassie import calibrate_threshold

                turd = settings.run_dir("lassie") / f"{settings.run_name}.turd"
                lassie_min = calibrate_threshold(turd, args.lassie_min_mad)
            print(f"Figure: {plot(settings, lassie_min=lassie_min, qseek_min=args.qseek_min_semblance, suffix=args.suffix, jitter=not args.no_jitter)}")
        case "summary":
            from lassie3.compare import summarise

            print(json.dumps(summarise(settings), indent=2))
        case "review":
            from lassie3.review import render, review

            print(render(review(settings)))
        case "sweep":
            from lassie3.reference import load as load_reference
            from lassie3.sweep import render as render_sweep, sweep

            rows = sweep(settings, args.mad, run=not args.no_run)
            print(render_sweep(rows, len(load_reference(settings))))
        case "ssst-terms":
            from lassie3.ssst import StationTerms, TermsSettings, derive
            from lassie3.ssst import render as render_terms

            terms_settings = TermsSettings(
                radius=args.radius,
                prior_weight=args.prior_weight,
                outlier_level=args.outlier_level,
                min_probability=args.min_probability,
            )
            path = derive(settings, terms_settings, days=args.train_days, limit=args.limit, force=args.force)
            print(render_terms(StationTerms.load(path)))
            held_out = path.parent / "held-out.json"
            if held_out.exists():
                print(f"Held-out check on {settings.day_str}: {held_out.read_text()}")
            print(f"Station terms: {path}")
        case "compare-ssst":
            from lassie3.ssst_compare import compare
            from lassie3.ssst_compare import render as render_ssst

            print(render_ssst(compare(settings)))

    return 0


if __name__ == "__main__":
    sys.exit(main())
