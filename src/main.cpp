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
vector<vector<double>> run_alg(string seq_a, string seq_b);
void print_matrix(vector<vector<double>> matrix, string file_name);
vector<string> read_files(vector<string> files);
map<string, string> parse_cl(int argc, char **argv);
void print_usage(char **argv);
void tally_diags(vector<vector<double>> matrix, string output_file);

//Mismatch score should be positive, gap penalty should be positive, as they are SUBTRACTED in the alg.
int mismatch_score = 1;
int gap_penalty = 0;

vector<vector<double>> run_alg(string seq_a, string seq_b)
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
  vector<vector<double>> score_matrix;
  vector<vector<int>> lower_traceback, upper_traceback;
  // Make all three matrices the correct size.
  score_matrix.resize(len_a + 1);
  lower_traceback.resize(len_a + 1);
  upper_traceback.resize(len_a + 1);
  for (auto &i : score_matrix)
  {
    i.resize(len_b + 1);
  }
  for (auto &i : lower_traceback)
  {
    i.resize(len_b + 1);
  }
  for (auto &i : upper_traceback)
  {
    i.resize(len_b + 1);
  }
  for (int i = 0; i <= len_a; i++)
  {
    for (int j = 0; j <= len_b; j++)
    {
      score_matrix[i][j] = 0.;
    }
  }
  // Fill maxtrices
  for (int i = 1; i <= len_a; i++)
  {
    for (int j = 1; j <= len_b; j++)
    {
      vector<double> values_to_compare;
      // If on the diagonal itll make the compare array [-1,-1,-1,0] such that 0 always wins saying this is the start of a new subsequence.
      if (i == j)
      {
        values_to_compare = {-1, -1, -1};
      }
      // If not on the diag do the regular smith waterman calculation.
      else
      {
        values_to_compare.push_back(score_matrix[i - 1][j - 1] + (seq_a[i - 1] == seq_b[j - 1] ? 1. : -mismatch_score));
        values_to_compare.push_back(score_matrix[i - 1][j] - gap_penalty);
        values_to_compare.push_back(score_matrix[i][j - 1] - gap_penalty);
      }
      // 0 is always in the running so it's not in a control statement
      values_to_compare.push_back(0.);
      // Find the max element and save it to the score matrix
      score_matrix[i][j] = *max_element(values_to_compare.begin(), values_to_compare.end());
      // Figure out which index is the max, therefore which score the max came from, to assist with the traceback... not importatn here except for completeness
      _int64 traceback_location = find(values_to_compare.begin(), values_to_compare.end(), H[i][j]) - values_to_compare.begin();
      switch (traceback_location)
      {
      case 0: // Match or mismatch
        lower_traceback[i][j] = i - 1;
        upper_traceback[i][j] = j - 1;
        break;
      case 1: // Deletion in A
        lower_traceback[i][j] = i - 1;
        upper_traceback[i][j] = j;
        break;
      case 2: // Deletion in B
        lower_traceback[i][j] = i;
        upper_traceback[i][j] = j - 1;
        break;
      case 3: // Start of new subseq
        lower_traceback[i][j] = i;
        upper_traceback[i][j] = j;
        break;
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

void tally_diags(vector<vector<double>> matrix, string output_file = "")
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
      output_file_stream << i << ": " << diag_sum << endl;
    }
  }
  output_file_stream.close();
}

void dump_matrix(vector<vector<double>> matrix, string output_file = "")
{
  // Loop over the matrix to find the largest value so the other values can be padded and give a nice output.
  // Start the max out at negative infinity so anything is larger.
  double max(-std::numeric_limits<double>::infinity());
  for (auto &i : matrix)
  {
    for (auto &j : i)
    {
      if (j > max)
      {
        max = j;
      }
    }
  }
  // Acutal dumping occurs here. Open the fille, write out the numbers in a sensical manner.
  string max_string = to_string(max);
  size_t pad_size = max_string.size();
  ofstream output_file_stream;
  output_file_stream.open(output_file, ofstream::out);
  if (output_file_stream.is_open())
  {
    for (auto &i : matrix)
    {
      for (auto &j : i)
      {
        size_t num_pad = pad_size - to_string(j).size();
        string pad = string(num_pad, ' ');
        output_file_stream << pad << j << " ";
      }
      output_file_stream << endl;
    }
  }
}

int main(int argc, char **argv)
{
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
    auto time = chrono::system_clock::now();
    auto int_time = chrono::system_clock::to_time_t(time);
    string string_time = ctime(&int_time);
    string_time = string_time.substr(0, string_time.find("\n"));
    string_time.replace(string_time.find("  "), 2, " ");
    replace(string_time.begin(), string_time.end(), ' ', '_');
    replace(string_time.begin(), string_time.end(), ':', '.');
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
  // Program 100% guarenteed to end up here.
  return 0;
}
