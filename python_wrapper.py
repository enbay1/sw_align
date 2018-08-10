# Short example and implementation of wrapping sw_align in python.
# Requires changes to run on linux, remove .exe in line 15.
import shlex
import subprocess


def main():
    file_1 = "./resources/test_short.fasta"
    file_2 = "./resources/test_short.fasta"
    output_file_name = "output.txt"
    matrix_dump_file_name = "dump_file.txt"
    gap_open = 25
    gap_extend = 1
    diagonal = "yes"
    command = "./build/sw_align.exe {} {} --output {} --dump {} --open {} --extend {} --diagonal {}".format(
        file_1, file_2, output_file_name, matrix_dump_file_name, gap_open, gap_extend, diagonal)
    print(subprocess.check_output(command).decode('utf-8'))
    return


if __name__ == "__main__":
    main()
