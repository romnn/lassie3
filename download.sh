#!/usr/bin/env bash
# Downloads StationXML (level=response) and miniSEED for the ICDP-EGER virtual
# network (_ICDPEGER), plus nearby open stations, for 2024-03-16 .. 2024-03-30
# (14 days).
# Restricted streams need an EIDA token: EIDA_TOKEN=~/.eidatoken ./download.sh
set -uo pipefail

START=2024-03-16
DAYS=14
ROOT="$(cd "$(dirname "$0")" && pwd)/data_icdp_eger"
GEOFON=https://geofon.gfz.de/fdsnws
BGR=https://eida.bgr.de/fdsnws

# node|net|stations
STATIONS=(
  "$GEOFON|WB|KOC KRC LBC NKC SKC STC VAC"
  "$GEOFON|6A|LWS00 LWSA1 LWSA2 LWSA3"
  "$GEOFON|1D|STC00"
  "$BGR|SX|ROHR TANN WERN GUNZ MULD WERD TRIB"
  # Open stations outside the virtual network, co-located with WB.NKC and 6A.LWS00.
  "$GEOFON|CZ|NKC"
  "$BGR|GQ|LNDWU"
)

mkdir -p "$ROOT/stationxml" "$ROOT/mseed" "$ROOT/logs"
end_day=$(date -j -v+${DAYS}d -f %Y-%m-%d "$START" +%Y-%m-%d)

# Fetches one URL to a file; keeps only HTTP 200 bodies.
fetch() {
  local url=$1 out=$2 code
  [[ -s $out ]] && { echo "skip $out"; return; }
  local opts=()
  if [[ $url == "$GEOFON/dataselect/"* && ${#AUTH[@]} -gt 0 ]]; then
    url=${url/\/query\?/\/queryauth?}; opts=("${AUTH[@]}")
  fi
  code=$(curl -sS --retry 5 --retry-delay 10 ${opts[@]+"${opts[@]}"} -o "$out.part" -w '%{http_code}' "$url")
  if [[ $code == 200 ]]; then
    mv "$out.part" "$out"; echo "ok   $code $(du -h "$out" | cut -f1) $out"
  else
    echo "fail $code $out: $(head -c 200 "$out.part" | tr '\n' ' ')"; rm -f "$out.part"
  fi
}

# With EIDA_TOKEN set, trade the token for temporary GEOFON credentials so
# restricted streams go through /queryauth with HTTP digest auth.
AUTH=()
if [[ -n ${EIDA_TOKEN:-} ]]; then
  creds=$(curl -sS --data-binary "@$EIDA_TOKEN" "$GEOFON/dataselect/1/auth") || exit 1
  [[ $creds == *:* ]] || { echo "token rejected: $creds"; exit 1; }
  AUTH=(--digest -u "$creds")
fi

for entry in "${STATIONS[@]}"; do
  IFS='|' read -r node net stas <<<"$entry"
  for sta in $stas; do
    fetch "$node/station/1/query?network=$net&station=$sta&level=response&starttime=$START&endtime=$end_day" \
      "$ROOT/stationxml/$net.$sta.xml"
  done
done

for entry in "${STATIONS[@]}"; do
  IFS='|' read -r node net stas <<<"$entry"
  for sta in $stas; do
    mkdir -p "$ROOT/mseed/$net/$sta"
    for ((i = 0; i < DAYS; i++)); do
      d0=$(date -j -v+${i}d -f %Y-%m-%d "$START" +%Y-%m-%d)
      d1=$(date -j -v+$((i + 1))d -f %Y-%m-%d "$START" +%Y-%m-%d)
      fetch "$node/dataselect/1/query?network=$net&station=$sta&starttime=${d0}T00:00:00&endtime=${d1}T00:00:00" \
        "$ROOT/mseed/$net/$sta/$net.$sta.$d0.mseed"
    done
  done
done
