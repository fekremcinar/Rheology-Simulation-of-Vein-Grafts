# Converting mov to mp4

```shell
ffmpeg -i "Screen Recording 2026-02-28 at 12.24.28.mov" \
  -c:v libx264 -preset slow -crf 26 \
  -c:a aac -b:a 96k \
  "heart_beat_in_vessel.mp4"
```

