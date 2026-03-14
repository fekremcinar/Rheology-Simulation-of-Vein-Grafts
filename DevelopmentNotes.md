## Converting mov to mp4

```shell
ffmpeg -i "Screen Recording 2026-02-28 at 12.24.28.mov" \
  -c:v libx264 -preset slow -crf 26 \
  -c:a aac -b:a 96k \
  "heart_beat_in_vessel.mp4"
```

## Reset ParaView

To reset ParaView settings on macOS, first quit ParaView, then run the following commands in the terminal:

```shell
rm -rf ~/Library/Preferences/ParaView
rm -rf ~/.config/ParaView
rm -rf ~/Library/Preferences/ParaView   # already done
rm -rf ~/Library/Application\ Support/ParaView
defaults delete org.paraview.ParaView 2>/dev/null || true
```

Reopen ParaView — the macros will be gone from the menu. Re-register with Tools → Macros → Add new macro.

## Run Experiments

```bash
docker run -it --rm \
  -u "$(id -u)":"$(id -g)" \
  -v "$HOME/Rheology-Simulation-of-Vein-Grafts":/work \
  opencfd/openfoam-default:2512
```

Inside the container:
```bash
cp -r /work/experiments/02_heartbeat_laminar /work/run/
cd /work/run/02_heartbeat_laminar
touch 02_heartbeat_laminar.foam
blockMesh
checkMesh
pimpleFoam
```


```bash
cp -r /work/experiments/03_vessel_junction /work/run/
cd /work/run/03_vessel_junction
touch 03_vessel_junction.foam
blockMesh
checkMesh
pimpleFoam
```