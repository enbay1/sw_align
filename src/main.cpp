#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <string>
#include <vector>
using namespace std;
// for color use \x1B[31m
// Forward declarations
struct score_matrix {
  vector<vector<int> > matrix;
  map<char, int> indices;
};
vector<vector<long long int> > run_alg(string seq_a, string seq_b,
                                       int open_gap_penalty,
                                       int extend_gap_penalty, int diagonal);
void print_matrix(vector<vector<long long int> > matrix, string file_name);
vector<string> read_files(vector<string> files);
map<string, string> parse_cl(int argc, char **argv);
void print_usage(char **argv);
void tally_diags(vector<vector<long long int> > matrix, string output_file);
score_matrix generate_matrix();
string to_lower(string to_lower);
// Converts a string to all lowercase. Mainly  used for CLI comparison to be
// case-insensitive.
string to_lower(string to_lower) {
  string lower = to_lower;
  for (int i = 0; i < lower.size(); i++)
    if (lower[i] >= 'A' && lower[i] <= 'Z') lower[i] = tolower(lower[i]);
  return lower;
}
// Runs the Smith-Waterman algorithm with seq_a vs seq_b with gap penalties.
// Diagonal determines whether to keep the major diagonal or not.
// If diagonal is set to 0 the main diagonal will be zeroed out.
vector<vector<long long int> > run_alg(string seq_a, string seq_b,
                                       int open_gap_penalty,
                                       int extend_gap_penalty, int diagonal) {
  // todo see if this is actually required.
  if (seq_b.length() < seq_a.length()) {
    cout << "Sequences have been swapped because seq a is longer than seq b. "
            "Consider passing files in opposite order."
         << endl;
    seq_a.swap(seq_b);
  }
  // get the actual lengths of the sequences
  size_t len_a = seq_a.length();
  size_t len_b = seq_b.length();
  // initialize matrix for scores.
  vector<vector<long long int> > score_matrix, seq_b_indel_matrix,
      seq_a_indel_matrix;
  // Make all three matrices the correct size.
  score_matrix.resize(len_a + 1);
  seq_b_indel_matrix.resize(len_a + 1);
  seq_a_indel_matrix.resize(len_a + 1);
  for (auto &i : score_matrix) {
    i.resize(len_b + 1);
  }
  for (auto &i : seq_b_indel_matrix) {
    i.resize(len_b + 1);
  }
  for (auto &i : seq_a_indel_matrix) {
    i.resize(len_b + 1);
  }
  // Fill maxtrices
  struct score_matrix dna_score_matrix = generate_matrix();
  for (int i = 1; i <= len_a; i++) {
    for (int j = 1; j <= len_b; j++) {
      // Grab characters to be compared
      char char_a = seq_a[i - 1];
      char char_b = seq_b[j - 1];
      // Get the chars' index in matrix
      int char_a_index = dna_score_matrix.indices[char_a];
      int char_b_index = dna_score_matrix.indices[char_b];
      // get score
      int score = dna_score_matrix.matrix[char_a_index][char_b_index];
      // set the next cell on the diagonal.
      seq_b_indel_matrix[i][j] =
          max((seq_b_indel_matrix[i - 1][j] - extend_gap_penalty),
              (score_matrix[i - 1][j] - open_gap_penalty));
      seq_a_indel_matrix[i][j] =
          max((seq_a_indel_matrix[i][j - 1] - extend_gap_penalty),
              (score_matrix[i][j - 1] - open_gap_penalty));
      score_matrix[i][j] = max(
          (long long)0,
          (max((score_matrix[i - 1][j - 1] + score),
               max((seq_b_indel_matrix[i][j]), (seq_a_indel_matrix[i][j])))));
      // This is where the diagonal flag is used. It's set to 0 after to retain
      // traceback info in the two indel matrices.
      if (i == j && !diagonal) {
        score_matrix[i][j] = 0;
      }
    }
  }
  return score_matrix;
}
// Read contents of two passed FASTA files
vector<string> read_files(vector<string> files) {
  vector<string> return_vector;
  string line;
  for (auto &file : files) {
    // Initiate file_contents which holds file contents
    vector<string> file_contents;
    // Try to open the file
    ifstream input_file(file, ios::in);
    // Check the file actually opened
    if (!input_file) {
      cerr << "Unable to open file " + file << endl;
      exit(1);
    } else {
      // throw DNA lines into a vector
      while (getline(input_file, line)) {
        if (line.find(">") == string::npos) {
          file_contents.push_back(line);
        }
      }
      input_file.close();
    }
    // Add all the fasta lines into one line
    return_vector.push_back(
        accumulate(file_contents.begin(), file_contents.end(), string("")));
  }
  return return_vector;
}
// Parse the command line.
map<string, string> parse_cl(int argc, char **argv) {
  vector<string> args;
  for (int i = 1; i < argc; i++) {
    args.push_back(argv[i]);
  }
  // See if the user asked for help.
  // no use running other things if use doesn't know what they're doing.
  for (const auto &i : args) {
    if (to_lower(i).compare("--help") == 0 || to_lower(i).compare("-h") == 0) {
      print_usage(argv);
    }
  }
  // 3 is <program name> <file 1> <file 2> and therefore minimum viable call. <
  // 3 is useless.
  if (argc < 3) {
    cerr << "No fasta files specified. Please pass 2 fasta files to align."
         << endl;
    exit(1);
  }
  map<string, string> parameters;
  // Optional paramaters need defaults.
  parameters.emplace("open_gap", "25");
  parameters.emplace("extend_gap", "1");
  parameters.emplace("diagonal", "1");
  parameters.emplace("output", "");
  parameters.emplace("dump", "");

  // Check that the files are of type .fasta
  // Files have to be in args[0] and [1] as they're the 1st and 2nd CLI params.
  for (string file : {args[0], args[1]}) {
    if (file.find(".fasta") != string::npos &&
        file.find(".FASTA") != string::npos) {
      cerr << "File " + file + " not of type .fasta" << endl;
      exit(1);
    }
  }
  // Store the two files in a map
  parameters.emplace("file_1", args[0]);
  parameters.emplace("file_2", args[1]);
  // Check for output names and dump option.
  if (args.size() > 12) {
    cerr << "Too many arguments specified." << endl;
    exit(1);
  }
  // This and the line that checks for less than three only exclude if there are
  // exactly three, which doesn't really require CLI parsing.
  if (args.size() > 3) {
    for (const auto &i : args) {
      if (i == args.back()) {
        break;
      }
      auto look_ahead = &i;
      look_ahead++;
      // Look for -o or --output for output file
      if (to_lower(i) == "-o" || to_lower(i) == "-output") {
        if (look_ahead[0] != "-") {
          parameters.find("output")->second = *look_ahead;
          continue;
        } else {
          cerr << "No output file name found after -o or --output flag. Using "
                  "generated one."
               << endl;
          parameters.emplace("output", "");
        }
        // Look for -d or --dump for matrix dump file
      } else if (to_lower(i) == "-d" || to_lower(i) == "--dump") {
        if (look_ahead[0] != "-") {
          parameters.find("dump")->second = *look_ahead;
          continue;
        } else {
          cerr << "No dump file name found after -d or --dump flag.  Using "
                  "generated one."
               << endl;
        }
        // Look for --open for gap open penalty
      } else if (to_lower(i) == "--open") {
        if (look_ahead[0] != "-") {
          parameters.find("open_gap")->second = *look_ahead;
          continue;
        } else {
          cerr << "No open gap penalty found after --open.  Using 25" << endl;
        }
        // Look for --extend for gap extension penalty
      } else if (to_lower(i) == "--extend" || to_lower(i) == "-e") {
        if (look_ahead[0] != "-") {
          parameters.find("extend_gap")->second = *look_ahead;
          continue;
        } else {
          cerr << "No extend gap penalty found after --open.  Using 1" << endl;
        }
        // Look for --diag or --diagonal for diagonal retention
      } else if (to_lower(i) == "--diag" || to_lower(i) == "--diagonal") {
        if (look_ahead[0] != "-") {
          if (to_lower(*look_ahead) == "no" || *look_ahead == "0") {
            parameters.find("diagonal")->second = "0";
          }
        } else {
          cerr << "No open gap penalty found after --open.  Using 25" << endl;
        }
      }
    }
  }
  return parameters;
}

void print_usage(char **argv) {
  cerr << "usage:"
       << " " << argv[0] << " "
       << "fasta_file_1"
       << " "
       << "fasta_file_2"
       << " "
       << "[-o output_file]"
       << " "
       << "[-d file_name]"
       << " "
       << "[--open gap_open_penalty]"
       << " "
       << "[-e gap_extend_penalty]"
       << " "
       << "[--diag keep_major_diagonal_intact]"
       << " "
       << "[-h]" << endl;
  cerr << endl;
  cerr << "Options:" << endl;
  cerr << "  fasta_file_1                 Name of fasta file 1 to be read in."
       << endl;
  cerr << "  fasta_file_2                 Name of fasta file 2 to be read in."
       << endl;
  cerr << "  -o, --output FILE            Name of output file. If one is not "
          "specified one will be generated."
       << endl;
  cerr << "  -d, --dump FILE              Name of file to which scoring matrix "
          "is written. If one is not specified one will be generated."
       << endl;
  cerr << "  --open                       Set the gap open penalty. A positive "
          "integer. [Default: 25]."
       << endl;
  cerr << "  -e, --extend (positive int)  Set the gap extend penalty. A "
          "positive integer. [Default: 1]."
       << endl;
  cerr << "  --diag, --diagonal [0|1]     Set whether to keep the major "
          "diagonal or not [Default: 1, keep diagonal]."
       << endl;
  cerr << "  -h, --help                   Show this message." << endl;
  exit(0);
}

void tally_diags(vector<vector<long long int> > matrix,
                 string output_file = "") {
  matrix.erase(matrix.begin());
  for (auto &i : matrix) {
    i.erase(i.begin());
  }
  ofstream output_file_stream;
  output_file_stream.open(output_file, ofstream::out);
  if (output_file_stream.is_open()) {
    // i goes across the top row of the matrix
    for (int i = 0; i < matrix[0].size(); i++) {
      long long int diag_sum = 0;
      // j is the running diagonal index  that starts the "column"
      int j = i;
      // k is the running diagonal index that goes from 0 to the end.
      int k = 0;
      while (j < matrix[0].size() && k < matrix.size()) {
        diag_sum += matrix[k][j];
        j++;
        k++;
      }
      output_file_stream << i << ":" << diag_sum * 2 << endl;
    }
    /*
    for (int i = 1; i < matrix.size(); i++) {
      double diag_sum = 0;
      // j is the running diagonal index  that starts the "column"
      int j = i;
      // k is the running diagonal index that goes from 0 to the end.
      int k = 0;
      while (j < matrix.size() && k < matrix[0].size()) {
        diag_sum += matrix[j][k];
        j++;
        k++;
      }
      output_file_stream << "-" << i << "\t" << diag_sum << endl;
    } */
  }
  output_file_stream.close();
}

void dump_matrix(vector<vector<long long int> > matrix,
                 string output_file = "") {
  // Loop over the matrix to find the largest value so the other values can be
  // padded and give a nice output. Start the max out at negative infinity so
  // anything is larger.
  int max = std::numeric_limits<int>::min();
  int min = std::numeric_limits<int>::max();
  for (auto &i : matrix) {
    for (auto &j : i) {
      if (j > max) {
        max = (int)j;
      }
      if (j < min) {
        min = int(j);
      }
    }
  }
  // If there's a neg number longer (including '-') than the longest pos number,
  // use the neg number to pad the output.
  auto pad_size = to_string(max).length();
  if (to_string(min).length() > pad_size) {
    pad_size = to_string(min).length();
  }
  // Actual dumping occurs here. Open the fille, write out the numbers in a
  // sensible manner.
  ofstream output_file_stream;
  output_file_stream.open(output_file, ofstream::out);
  if (output_file_stream.is_open()) {
    for (int i = 1; i < matrix.size(); i++) {
      for (int j = 1; j < matrix[i].size(); j++) {
        auto j_string = to_string(matrix[i][j]);
        size_t amount_to_pad = pad_size - j_string.length();
        string pad = string(amount_to_pad, ' ');
        output_file_stream << pad << (int)matrix[i][j] << " ";
      }
      output_file_stream << endl;
    }
  }
}

score_matrix generate_matrix() {
  score_matrix dna_score_matrix;
  dna_score_matrix.indices['A'] = 0;
  dna_score_matrix.indices['T'] = 1;
  dna_score_matrix.indices['G'] = 2;
  dna_score_matrix.indices['C'] = 3;
  dna_score_matrix.matrix = {
      {5, -4, -4, -4}, {-4, 5, -4, -4}, {-4, -4, 5, -4}, {-4, -4, -4, 5}};
  return dna_score_matrix;
}

int main(int argc, char **argv) {
  chrono::high_resolution_clock::time_point t1 =
      chrono::high_resolution_clock::now();
  // Parse CLI
  map<string, string> cl_params = parse_cl(argc, argv);
  // Read files from CLI
  vector<string> file_contents =
      read_files({cl_params["file_1"], cl_params["file_2"]});
  // Run the algorithm
  auto matrix =
      run_alg(file_contents[0], file_contents[1], stoi(cl_params["open_gap"]),
              stoi(cl_params["extend_gap"]), stoi(cl_params["diagonal"]));
  // Check to see if output and dump files were specified. If they weren't
  // assign them regardless of use. This way if they're not specified they're
  // the same name and will sort nicely in the OS.
  if (cl_params["output"] == "" || cl_params["dump"] == "") {
    time_t t = chrono::system_clock::to_time_t(chrono::system_clock::now());
    char buf[100] = {0};
    strftime(buf, sizeof(buf), "%Y.%m.%d_%H.%M.%S", localtime(&t));
    auto string_time = string(buf);
    // But don't overwrite someones' hard-input name! That's rude.
    if (cl_params["output"] == "") {
      cl_params["output"] = "sw_align_" + string_time + "_output.txt";
    }
    if (cl_params["dump"] == "") {
      cl_params["dump"] = "sw_align_" + string_time + "_matrix.txt";
    }
  }
  // Sum the diagonals, the entire point of this program.
  tally_diags(matrix, cl_params["output"]);
  // "helpful" messages to end users!
  cout << "Results saved to " << cl_params["output"] << endl;
  // See if the end user put "-d" anywhere, if they did dump ze matrix!
  for (int i = 1; i < argc; i++) {
    if (string(argv[i]).compare("-d") == 0 ||
        string(argv[i]).compare("-dump") == 0) {
      dump_matrix(matrix, cl_params["dump"]);
      // Always tell them where you put their file!
      cout << "Matrix saved to " << cl_params["dump"] << endl;
    }
  }
  chrono::high_resolution_clock::time_point t2 =
      chrono::high_resolution_clock::now();
  auto duration = chrono::duration_cast<chrono::milliseconds>(t2 - t1).count();
  cout << "SW align completed in " << duration << "ms" << endl;
  // Program 100% guaranteed to end up here.
  return 0;
}
