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