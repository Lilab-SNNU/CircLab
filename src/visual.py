import sys
import os
import argparse
from tabulate import tabulate


def generate_help_message():
    """Generate help message with a table of functions and descriptions"""
    table = [
        ['ClassByType', 'Visualizes circRNA classification by type. -i <input_file> -o <output_file>'],
        ['showJunction', 'Displays junction reads for circRNA. -i <input_file> -n <number> -o <output_file>'],
        ['showTissueDistribution',
         'Shows distribution of circRNA across tissues. -i <input_file> [contain tissue and circRNA_type column] -o <output_file>'],
        ['showLength', 'Displays length distribution of circRNA. -i <input_file> -l <length> -o <output_file>'],
        ['showLength_barplot', 'Creates a bar plot for circRNA lengths. -i <input_file> -l <length> -o <output_file>'],
        ['showCodon',
         'can help to display triplet codon composition frequency. -i <input_file> -n <number> -o <output_file>']
    ]

    help_message = """
Visual - Choose a visualization function to execute

Available Functions:
"""
    help_message += tabulate(table, headers=['Function Name', 'Description'], tablefmt='grid')
    return help_message


def main(args=None):
    parser = argparse.ArgumentParser(
        description=generate_help_message(),
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument('-f', '--function',
                        required=True,
                        help='Specify the function to execute (e.g., showTissueDistribution)')

    parser.add_argument('-i', '--input',
                        help='Input file for the function (required for showLength and showTissueDistribution)')

    parser.add_argument('-l', '--length',
                        type=int,
                        help='Length parameter for the function (required for showLength and showLength_barplot)')

    parser.add_argument('-o', '--output',
                        help='Output file for the function results (e.g., graphs)')

    parser.add_argument('-n', '--max_splice',
                        type=int,
                        default=20,
                        help='Maximum number of splice-signal types to display in barplot (default 20)')

    parsed_args = parser.parse_args(args)
    function_map = {
        "ClassByType": "visual/ClassByType.R",
        "showJunction": "visual/showJunction.R",
        "showTissueDistribution": "visual/showTissueDistribution.R",
        "showLength": "visual/showLength.R",
        "showLength_barplot": "visual/showLength_barplot.R",
        "showCodon": "visual/showCodon.R"
    }


    r_script = function_map.get(parsed_args.function)
    if not r_script:
        print(f"Error: Unknown function {parsed_args.function}.\n")
        print(generate_help_message())
        sys.exit(2)

    if parsed_args.function == "showCodon":
        if not all([parsed_args.input, parsed_args.output]):
            print("Error: -i <input_file> and -o <output_file> are required for showCodon.\n")
            sys.exit(2)

        r_command = f"Rscript {r_script} {parsed_args.input} {parsed_args.max_splice} {parsed_args.output}"

    elif parsed_args.function == "showJunction":
        if not all([parsed_args.input, parsed_args.output]):
            print("Error: -i <input_file> and -o <output_file> are required for showJunction.\n")
            sys.exit(2)

        r_command = f"Rscript {r_script} {parsed_args.input} {parsed_args.output} {parsed_args.max_splice}"

    elif parsed_args.function == "ClassByType":
        if not all([parsed_args.input, parsed_args.output]):
            print("Error: -i <input_file> and -o <output_file> are required for ClassByType.\n")
            sys.exit(2)

        r_command = f"Rscript {r_script} {parsed_args.input} {parsed_args.output}"

    elif parsed_args.function == "showLength":
        if not all([parsed_args.input, parsed_args.length, parsed_args.output]):
            print("Error: -i <input_file>, -l <length>, and -o <output_file> are required for showLength.\n")
            sys.exit(2)

        r_command = f"Rscript {r_script} {parsed_args.input} {parsed_args.length} {parsed_args.output}"

    elif parsed_args.function == "showLength_barplot":
        if not all([parsed_args.input, parsed_args.length, parsed_args.output]):
            print("Error: -i <input_file>, -l <length>, and -o <output_file> are required for showLength_barplot.\n")
            sys.exit(2)

        r_command = f"Rscript {r_script} {parsed_args.input} {parsed_args.length} {parsed_args.output}"

    else:
        if parsed_args.params:
            r_command = f"Rscript {r_script} {parsed_args.params}"
        else:
            print(f"Error: No additional parameters provided for {parsed_args.function}.")
            sys.exit(2)

    print(f"Executing: {r_command}")

    try:
        result = os.system(r_command)
        if result != 0:
            print(f"Error executing R script. Exit code: {result}")
            sys.exit(2)
    except Exception as e:
        print(f"An error occurred: {e}")
        sys.exit(2)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(generate_help_message())
        sys.exit(1)
    main()
