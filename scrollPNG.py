#!/usr/bin/env python3
import argparse
from math import ceil
from PIL import Image

def make_scrolling_gif(
    input_path,
    output_path,
    window_width=500,
    fps=30,
    duration=20.0,
):
    """
    Make a scrolling window animation over a wide image.

    - input_path:  source PNG (e.g. 100 x 10000)
    - output_path: output GIF
    - window_width: width of the visible window in pixels
    - fps: frames per second
    - duration: approximate total time in seconds
    """
    img = Image.open(input_path).convert("RGB")
    full_w, full_h = img.size

    if window_width > full_w:
        print(f"Window width {window_width} > image width {full_w}, using {full_w}")
        window_width = full_w

    total_dx = full_w - window_width

    # No scrolling possible: just save one frame
    if total_dx <= 0:
        frames = [img]
        ms_per_frame = int(1000 / fps)
    else:
        # Target number of frames for desired duration
        target_frames = max(2, int(fps * duration))

        # Step size in pixels per frame so we roughly cover total_dx in target_frames
        # We want positions from x=0 to x=total_dx (inclusive)
        step = max(1, total_dx // (target_frames - 1))

        positions = list(range(0, total_dx + 1, step))
        frames = [
            img.crop((x, 0, x + window_width, full_h))
            for x in positions
        ]

        actual_frames = len(frames)
        # Adjust duration per frame to approach target duration
        ms_per_frame = int(duration * 1000 / actual_frames)

    # Save GIF
    frames[0].save(
        output_path,
        save_all=True,
        append_images=frames[1:],
        duration=ms_per_frame,
        loop=0,
        optimize=False,
    )
    print(
        f"Saved scrolling GIF to {output_path} "
        f"({len(frames)} frames, ~{ms_per_frame} ms/frame)"
    )


def main():
    parser = argparse.ArgumentParser(
        description="Create a horizontal scrolling GIF from a wide PNG."
    )
    parser.add_argument("input", help="Input PNG path")
    parser.add_argument("output", help="Output GIF path")
    parser.add_argument(
        "--window-width",
        type=int,
        default=500,
        help="Width of the scrolling window in pixels (default: 500)",
    )
    parser.add_argument(
        "--fps",
        type=int,
        default=30,
        help="Frames per second (default: 30)",
    )
    parser.add_argument(
        "--duration",
        type=float,
        default=20.0,
        help="Approximate total duration in seconds (default: 20)",
    )
    args = parser.parse_args()

    make_scrolling_gif(
        args.input,
        args.output,
        window_width=args.window_width,
        fps=args.fps,
        duration=args.duration,
    )


if __name__ == "__main__":
    main()
