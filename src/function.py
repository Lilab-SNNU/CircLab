import sys


def print_help():
    help_message = """
    Usage: python function.py <subcommand> [options]

    Subcommands:
        miRNA           Analyze miRNA interactions.
        conservation    Analyze circRNA conservation across species.

    Examples:
        circlab function -miRNA -i input_file -o output_file
        circlab function -conservation -s species.txt -i input_folder -o output.txt -t 25 -f 25

    Use "circlab function <subcommand> --help" for more details on each subcommand.
    """
    print(help_message)


def main(argv=None):
    if not argv or len(argv) < 1:
        print("Error: Missing subcommand. Use 'miRNA' or 'conservation'.")
        print_help()
        return

    subcommand = argv[0]  

    if subcommand == "miRNA":
        try:
            from src import miRNA  
            miRNA.main(argv[1:])  
        except ImportError as e:
            print(f"Error: Failed to import 'miRNA' module. {e}")
    elif subcommand == "conservation":
        try:
            from src import conservation  
            conservation.main(argv[1:])
        except ImportError as e:
            print(f"Error: Failed to import 'conservation' module. {e}")
    elif subcommand in ["help", "--help", "-h"]:
        print_help()
    else:
        print(f"Error: Unknown subcommand '{subcommand}'. Valid options: 'miRNA', 'conservation'.")
        print_help()

if __name__ == "__main__":
    main(sys.argv[1:]) 
