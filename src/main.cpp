#include <iostream>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <filesystem>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <numeric>
#include <chrono>
using namespace std;

// Forward declarations.
struct score_matrix
{
  vector<vector<int>> matrix;
  map<char, int> indices;
};
vector<vector<int>> run_alg(string seq_a, string seq_b);
void print_matrix(vector<vector<int>> matrix, string file_name);
vector<string> read_files(vector<string> files);
map<string, string> parse_cl(int argc, char **argv);
void print_usage(char **argv);
void tally_diags(vector<vector<int>> matrix, string output_file);
score_matrix generate_matrix();

//Scores should be positive, as they are SUBTRACTED in the alg.
int open_gap_penalty = 25;
int extend_gap_penalty = 1;

vector<vector<int>> run_alg(string seq_a, string seq_b)
{

  if (seq_b.length() < seq_a.length())
  {
    cout << "Sequences have been swapped becuase seq a is longer than seq b. Consider passing files in opposite order." << endl;
    seq_a.swap(seq_b);
  }
  // get the actual lengths of the sequences
  size_t len_a = seq_a.length();
  size_t len_b = seq_b.length();
  // initialize matrix for scores, and traceback, I guess.
  vector<vector<int>> score_matrix;
  vector<vector<int>> IC, IR;
  // Make all three matrices the correct size.
  score_matrix.resize(len_a + 1);
  IC.resize(len_a + 1);
  IR.resize(len_a + 1);
  for (auto &i : score_matrix)
  {
    i.resize(len_b + 1);
  }
  for (auto &i : IC)
  {
    i.resize(len_b + 1);
  }
  for (auto &i : IR)
  {
    i.resize(len_b + 1);
  }
  // Fill maxtrices
  auto dna_score_matrix = generate_matrix();
  for (int row = 0; row < len_a; row++)
  {
    for (int column = 0; column < len_b; column++)
    {
      auto char_a = seq_a[row];
      auto char_b = seq_b[column];
      auto char_a_index = dna_score_matrix.indices[char_a];
      auto char_b_index = dna_score_matrix.indices[char_b];
      int score = dna_score_matrix.matrix[char_a_index][char_b_index];
      if (row == column)
      {
        score_matrix[row + 1][column + 1] = 0;
        IR[row + 1][column + 1] = max(score_matrix[row][column + 1] - open_gap_penalty, IR[row][column + 1] - extend_gap_penalty);
        IC[row + 1][column + 1] = max(score_matrix[row + 1][column] - open_gap_penalty, IC[row + 1][column] - extend_gap_penalty);
      }
      else
      {
        score_matrix[row + 1][column + 1] = score + max(0, max(score_matrix[row][column], max(IR[row][column], IC[row][column])));
        IR[row + 1][column + 1] = max(score_matrix[row][column + 1] - open_gap_penalty, IR[row][column + 1] - extend_gap_penalty);
        IC[row + 1][column + 1] = max(score_matrix[row + 1][column] - open_gap_penalty, IC[row + 1][column] - extend_gap_penalty);
      }
    }
  }
  return score_matrix;
}

vector<string> read_files(vector<string> files)
{
  vector<string> return_vector;
  string line;
  for (auto &file : files)
  {
    // Initiate file_contents which holds file contents
    vector<string> file_contents;
    // Try to open the file
    ifstream input_file(file, ios::in);
    // Check the file actually opened
    if (!input_file)
    {
      cerr << "Unable to open file " + file << endl;
      exit(1);
    }
    else
    {
      while (getline(input_file, line))
      {
        if (line.find(">") == string::npos)
        {
          file_contents.push_back(line);
        }
      }
      input_file.close();
    }
    return_vector.push_back(accumulate(file_contents.begin(), file_contents.end(), string("")));
  }
  return return_vector;
}

map<string, string> parse_cl(int argc, char **argv)
{
  if (argc < 3)
  {
    print_usage(argv);
  }
  vector<string> args;
  for (int i = 1; i < argc; i++)
  {
    args.push_back(argv[i]);
  }
  map<string, string> parameters;
  // See if the user asked for help.
  for (const auto &i : args)
  {
    if (i.compare("--help") == 0 || i.compare("-h") == 0)
    {
      print_usage(argv);
    }
  }
  // Check that the files are of type .fasta
  for (string file : {args[0], args[1]})
  {
    if (file.find(".fasta") != string::npos && file.find(".FASTA") != string::npos)
    {
      cerr << "File " + file + " not of type .fasta" << endl;
      ;
      print_usage(argv);
    }
  }
  // Store the two files in a map
  parameters.emplace("file_1", args[0]);
  parameters.emplace("file_2", args[1]);
  // Check for output names and dump option.
  if (args.size() > 7)
  {
    cerr << "Too many arguments specified." << endl;
    print_usage(argv);
  }
  if (args.size() > 3)
  {
    for (const auto &i : args)
    {
      if (i == args.back())
      {
        break;
      }
      auto look_ahead = &i;
      look_ahead++;
      if (i == "-o" || i == "-output")
      {
        if (look_ahead->find("-") != 0 && look_ahead->find("--") != 0)
        {
          parameters.emplace("output", *look_ahead);
          continue;
        }
        else
        {
          cerr << "No output file name found after -o or --output flag. Using generated one." << endl;
          parameters.emplace("output", "");
        }
      }
      if (i == "-d" || i == "-dump")
      {
        if (look_ahead->find("-") != 0 && look_ahead->find("--") != 0)
        {
          parameters.emplace("dump", *look_ahead);
          continue;
        }
        else
        {
          cerr << "No dump file name found after -d or --dump flag.  Using generated one." << endl;
          parameters.emplace("dump", "");
        }
      }
    }
  }
  return parameters;
}

void print_usage(char **argv)
{
  cout << "Entering print_usage." << endl;
  cerr << "Usage:"
       << " " << argv[0] << " "
       << "fasta_file_1"
       << " "
       << "fasta_file_2"
       << " "
       << "[-o output_file]"
       << " "
       << "[-d file_name]"
       << "[-h]" << endl;
  cerr << endl;
  cerr << "Options:" << endl;
  cerr << "  fasta_file_1      Name of fasta file 1 to be read in." << endl;
  cerr << "  fasta_file_1      Name of fasta file 2 to be read in." << endl;
  cerr << "  -h --help         Show this message." << endl;
  cerr << "  -o --output FILE  Name of output file. If one is not specified one will be generated." << endl;
  cerr << "  -d --dump FILE    Name of file to which scoring matrix is written. If one is not specified one will be generated." << endl;
  exit(0);
}

void tally_diags(vector<vector<int>> matrix, string output_file = "")
{
  ofstream output_file_stream;
  output_file_stream.open(output_file, ofstream::out);
  if (output_file_stream.is_open())
  {
    // i goes across the top row of the matrix
    for (int i = 0; i < matrix[0].size(); i++)
    {
      double diag_sum = 0;
      // j is the running diagonal index  that starts the "column"
      int j = i;
      // k is the running diagonal index that goes from 0 to the end.
      int k = 0;
      while (j < matrix[0].size() && k < matrix.size())
      {
        diag_sum += matrix[k][j];
        j++;
        k++;
      }
      output_file_stream << i << "\t" << diag_sum << endl;
    }
  }
  output_file_stream.close();
}

void dump_matrix(vector<vector<int>> matrix, string output_file = "")
{
  // Loop over the matrix to find the largest value so the other values can be padded and give a nice output.
  // Start the max out at negative infinity so anything is larger.
  int max = std::numeric_limits<int>::min();
  int min = std::numeric_limits<int>::max();
  for (auto &i : matrix)
  {
    for (auto &j : i)
    {
      if (j > max)
      {
        max = (int)j;
      }
      if (j < min)
      {
        min = int(j);
      }
    }
  }
  // If there's a neg number longer (including '-') than the longest pos number, use the neg number to pad the output.
  auto pad_size = to_string(max).length();
  if (to_string(min).length() > pad_size)
  {
    pad_size = to_string(min).length();
  }
  // Acutal dumping occurs here. Open the fille, write out the numbers in a sensical manner.
  ofstream output_file_stream;
  output_file_stream.open(output_file, ofstream::out);
  if (output_file_stream.is_open())
  {
    for (auto &i : matrix)
    {
      for (auto &j : i)
      {
        auto j_string = to_string(j);
        size_t amount_to_pad = pad_size - j_string.length();
        string pad = string(amount_to_pad, ' ');
        output_file_stream << pad << (int)j << " ";
      }
      output_file_stream << endl;
    }
  }
}

score_matrix generate_matrix()
{
  score_matrix dna_score_matrix;
  dna_score_matrix.indices['A'] = 0;
  dna_score_matrix.indices['T'] = 1;
  dna_score_matrix.indices['G'] = 2;
  dna_score_matrix.indices['C'] = 3;
  dna_score_matrix.matrix = {{5, -4, -4, -4},
                             {-4, 5, -4, -4},
                             {-4, -4, 5, -4},
                             {-4, -4, -4, 5}};
  return dna_score_matrix;
}

int main(int argc, char **argv)
{
  chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
  // Parse CLI
  map<string, string> cl_params = parse_cl(argc, argv);
  // Read files from CLI
  vector<string> file_contents = read_files({cl_params["file_1"], cl_params["file_2"]});
  // Run the algorithm
  auto matrix = run_alg(file_contents[0], file_contents[1]);
  // Check to see if output and dump files were specified. If they weren't assign them regardless of use.
  // This way if they're not specified they're the same name and will sort nicely in the OS.
  if (cl_params["output"] == "" || cl_params["dump"] == "")
  {
    time_t t = chrono::system_clock::to_time_t(chrono::system_clock::now());
    char buf[100] = {0};
    strftime(buf, sizeof(buf), "%Y.%m.%d_%H.%M.%S", localtime(&t));
    auto string_time = string(buf);
    // But don't overwrite someones' hard-input name! That's rude.
    if (cl_params["output"] == "")
    {
      cl_params["output"] = "sw_align_" + string_time + "_output.txt";
    }
    if (cl_params["dump"] == "")
    {
      cl_params["dump"] = "sw_align_" + string_time + "_matrix.txt";
    }
  }
  // Sum the diagonals, the entire point of this program.
  tally_diags(matrix, cl_params["output"]);
  // "helpful" messages to end users!
  cout << "Results saved to " << cl_params["output"] << endl;
  // See if the end ueser put "-d" anywhere, if they did dump ze matrix!
  for (int i = 1; i < argc; i++)
  {
    if (string(argv[i]).compare("-d") == 0 || string(argv[i]).compare("-dump") == 0)
    {
      dump_matrix(matrix, cl_params["dump"]);
      // Always tell them where you put their file!
      cout << "Matrix saved to " << cl_params["dump"] << endl;
    }
  }
  chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
  auto duration = chrono::duration_cast<chrono::milliseconds>(t2 - t1).count();
  cout << "SW align completed in " << duration << "ms" << endl;
  // Program 100% guarenteed to end up here.
  return 0;
}
