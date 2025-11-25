#!/usr/bin/env python3
import argparse
import numpy as np
from PIL import Image

def squeeze_black_columns(input_path, output_path, threshold=0):
    """
    Load an RGB PNG, remove columns that are entirely 'black',
    and save the squeezed image.

    threshold: treat any pixel with all channels <= threshold as black.
               default 0 = strictly (0,0,0)
    """
    img = Image.open(input_path).convert("RGB")
    arr = np.array(img)  # shape (H, W, 3)

    # A pixel is "black" if all channels <= threshold
    black_pixels = (arr <= threshold).all(axis=2)      # shape (H, W)
    black_cols    = black_pixels.all(axis=0)           # shape (W,)
    keep_cols     = ~black_cols                        # True for informative columns

    if not keep_cols.any():
        raise ValueError("All columns are black; nothing to keep.")

    squeezed = arr[:, keep_cols, :]                    # shape (H, W', 3)
    out_img = Image.fromarray(squeezed, mode="RGB")
    out_img.save(output_path)
    print(f"Squeezed image saved to {output_path} (width {squeezed.shape[1]})")


def main():
    parser = argparse.ArgumentParser(
        description="Remove all-black columns from an RGB PNG."
    )
    parser.add_argument("input", help="Input PNG path")
    parser.add_argument("output", help="Output PNG path")
    parser.add_argument(
        "--threshold",
        type=int,
        default=0,
        help="Max channel value to still treat as black (default: 0)"
    )
    args = parser.parse_args()
    squeeze_black_columns(args.input, args.output, args.threshold)


if __name__ == "__main__":
    main()
