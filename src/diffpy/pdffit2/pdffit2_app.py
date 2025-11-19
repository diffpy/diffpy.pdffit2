import argparse

from diffpy.pdffit2.version import __version__  # noqa


def main():
    parser = argparse.ArgumentParser(
        prog="diffpy.pdffit2",
        description=(
            "PDFfit2 - real space structure refinement program.\n\n"
            "For more information, visit: "
            "https://github.com/diffpy/diffpy.pdffit2/"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "--version",
        action="store_true",
        help="Show the program's version number and exit",
    )

    args = parser.parse_args()

    if args.version:
        print(f"diffpy.pdffit2 {__version__}")
    else:
        # Default behavior when no arguments are given
        parser.print_help()


if __name__ == "__main__":
    main()
