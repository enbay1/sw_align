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



vector< vector<double> > run_alg(string seq_a, string seq_b);
void print_matrix(vector< vector<double> > matrix, string file_name);
vector<string> read_files(vector<string> files);
map<string, string> parse_cl(int argc, char **argv);
void print_usage(char **argv);
void tally_diags(vector< vector<double> > matrix, string output_file);

//Mismatch score should be positive, gap penalty should be positive, as they are SUBTRACTED in the alg.
int mismatch_score = 1;
int gap_penalty = 0;


vector< vector<double> > run_alg(string seq_a, string seq_b){
  // get the actual lengths of the sequences
  size_t len_a = seq_a.length(); 
  size_t len_b = seq_b.length();
  // initialize matrix for scores, and traceback, I guess.
  vector< vector<double> > H;
  vector< vector<int> > I_i, I_j;
  // Make all three matrices the correct size.
  H.resize(len_a + 1);
  I_i.resize(len_a + 1);
  I_j.resize(len_a + 1);
  for (auto &i : H){
    i.resize(len_b + 1);
  }
  for (auto &i : I_i){
    i.resize(len_b + 1);
  }
  for (auto &i : I_j){
    i.resize(len_b + 1);
  }
  for (int i = 0; i <= len_a; i++){
    for (int j = 0; j <= len_b; j++)
    {
      H[i][j] = 0.;
    }
  }
  // Fill maxtrices
  for (int i = 1; i <= len_a; i++){
    for (int j = 1; j <= len_b; j++){
      vector<double> values_to_compare;
      if(i==j){
        for(i=0; i < 3; i++){
          values_to_compare.push_back(-1);
        }
      }
      else{
        values_to_compare.push_back(H[i - 1][j - 1] + (seq_a[i - 1] == seq_b[j - 1] ? 1. : -mismatch_score));
        values_to_compare.push_back(H[i - 1][j] - gap_penalty);
        values_to_compare.push_back(H[i][j - 1] - gap_penalty);
      }
      values_to_compare.push_back(0.);
      H[i][j] = *max_element(values_to_compare.begin(), values_to_compare.end());
      _int64 ind = find(values_to_compare.begin(), values_to_compare.end(), H[i][j]) - values_to_compare.begin();
      switch (ind){
      case 0: // Match or mismatch
        I_i[i][j] = i - 1;
        I_j[i][j] = j - 1;
        break;
      case 1: // Deletion in A
        I_i[i][j] = i - 1;
        I_j[i][j] = j;
        break;
      case 2: // Deletion in B
        I_i[i][j] = i;
        I_j[i][j] = j - 1;
        break;
      case 3: // 0
        I_i[i][j] = i;
        I_j[i][j] = j;
        break;
      }
    }
  }
  return H;
}

void print_matrix(vector< vector<double> > matrix, string file_name){
  for (int i = 1; i <= matrix.size(); i++){
    for (int j = 1; j <=  matrix[0].size(); j++){
      cout << matrix[i][j] << " ";
    }
    cout << endl;
  }
}

vector<string> read_files(vector<string> files){
    vector<string> return_vector;
    string line;
    for (auto &file: files){
        // Initiate file_contents which holds file dontents
        vector<string> file_contents;
        // Try to open the file
        ifstream input_file(file, ios::in);
        // Check the file actually opened
        if (!input_file) {
            cerr << "Unable to open file " + file << endl;
            exit(1);
        }
        else{
            while (getline(input_file, line)){
              if(line.find(">") == string::npos){
                file_contents.push_back(line);
              }
            }
            input_file.close();
        }
        //string contents = string(accumismatch_scorelate(file_contents.begin(), file_contents.end(), string(""))));
        return_vector.push_back(accumulate(file_contents.begin(), file_contents.end(), string("")));
    }
    return return_vector;
}

map<string, string> parse_cl(int argc, char **argv){
  vector<string> args;
  for(int i=1; i < argc; i++){
    args.push_back(argv[i]);
  }
    map<string, string> parameters;
    // See if the user asked for help.
    for(const auto &i: args){
      if(i.compare("--help") == 0 || i.compare("-h") == 0){
        print_usage(argv);
      }
    }
    // Check that the files are of type .fasta
    for (string file: {args[0], args[1]}){
        if (file.find(".fasta") != string::npos && file.find(".FASTA") != string::npos){
            cerr << "File "  + file + " not of type .fasta" << endl;;
            print_usage(argv);
        }
    }
    // Store the two files in a map
    parameters.emplace("file_1", args[0]);
    parameters.emplace("file_2", args[1]);
    // Check for output names and dump option.
    if (args.size() > 7){
      cerr << "Too many arguments specified." << endl;
      print_usage(argv);
    }
    if (args.size() > 3){
      for(const auto &i: args){
        if(i == args.back()){
          break;
        }
        auto look_ahead = &i;
        look_ahead++;
          if (i == "-o" || i == "-output"){
            if(look_ahead->find("-") != 0 && look_ahead->find("--") != 0){
              parameters.emplace("output", *look_ahead);
              continue;
            }
            else{
              cerr << "No output file name found after -o or --output flag. Using generated one." << endl;
              parameters.emplace("output", "");
            }
          }
          if (i == "-d" || i == "-dump"){
            if(look_ahead->find("-") != 0 && look_ahead->find("--") != 0){
              parameters.emplace("dump", *look_ahead);
              continue;
            }
            else{
              cerr << "No dump file name found after -d or --dump flag.  Using generated one." << endl;
              parameters.emplace("dump", "");
            }
          }
      }
    }
    return parameters;
}

void print_usage(char **argv){
  cerr << "Usage:" << " " << argv[0] << " " << "fasta_file_1" << " " << "fasta_file_2" << " " << "[-o output_file]" << " " << "[-d file_name]" << "[-h]" << endl;
  cerr << endl;
  cerr << "Options:" << endl;
  cerr << "  fasta_file_1      Name of fasta file 1 to be read in." << endl;
  cerr << "  fasta_file_1      Name of fasta file 2 to be read in."  << endl;
  cerr << "  -h --help         Show this message." << endl;
  cerr << "  -o --output FILE  Name of output file. If one is not specified one will be generated." << endl;
  cerr << "  -d --dump FILE    Name of file to which scoring matrix is written. If one is not specified one will be generated." << endl;
  exit(0);
}

void tally_diags(vector< vector<double> > matrix, string output_file = ""){
  if(output_file == ""){
    auto time = chrono::system_clock::now();
    auto int_time = chrono::system_clock::to_time_t(time);
    string string_time = ctime(&int_time);
    cout << endl;
    string_time = string_time.substr(0, string_time.find("\n"));
    string_time.replace(string_time.find("  "),2, " ");
    replace(string_time.begin(), string_time.end(), ' ', '_');
    replace(string_time.begin(), string_time.end(), ':', '.');
    output_file = "sw_align_" + string_time + ".txt";
  }
  ofstream output_file_stream;
  output_file_stream.open(output_file, ofstream::out);
  if(output_file_stream.is_open()){
    for(int i=0; i < matrix[0].size(); i++){
      double diag_sum = 0;
      auto j = i;
      auto k = 0;
      while(j != matrix[0].size()){
        diag_sum += matrix[k][j];
        j++;
        k++;
      }
      output_file_stream << i << ": " << diag_sum << endl;
    }
  }
  output_file_stream.close();
}

void dump_matrix(vector< vector<double> > matrix, string output_file = ""){
  matrix = {{1,2,4}};
  matrix.push_back({1,2,100});
  double max(-std::numeric_limits<double>::infinity());
  for(auto &i: matrix){
    for (auto &j: i){
      if(j > max){
        max = j;
      }
    }
  }
  string max_string = to_string(max);
  size_t pad_size = max_string.size();
  ofstream output_file_stream;
  output_file_stream.open(output_file, ofstream::out);
  if(output_file_stream.is_open()){
    for(auto &i: matrix){
      for (auto &j: i){
        size_t num_pad = pad_size - to_string(j).size();
        string pad = string(num_pad, ' ');
        output_file_stream << pad << j << " ";
      }
      output_file_stream << endl;
    }
  }
}

int main(int argc, char **argv){
    map<string, string> cl_params = parse_cl(argc, argv);
    vector<string> file_contents = read_files({cl_params["file_1"], cl_params["file_2"]});
    auto matrix = run_alg(file_contents[0], file_contents[1]);
    if(cl_params["output"] == "" || cl_params["dump"] == ""){
      auto time = chrono::system_clock::now();
      auto int_time = chrono::system_clock::to_time_t(time);
      string string_time = ctime(&int_time);
      string_time = string_time.substr(0, string_time.find("\n"));
      string_time.replace(string_time.find("  "),2, " ");
      replace(string_time.begin(), string_time.end(), ' ', '_');
      replace(string_time.begin(), string_time.end(), ':', '.');
      if(cl_params["output"] == ""){
        cl_params["output"] = "sw_align_" + string_time + "_output.txt";
      }
      if(cl_params["dump"] == ""){
        cl_params["dump"] = "sw_align_" + string_time + "_matrix.txt";
      }
    }
    tally_diags(matrix, cl_params["output"]);
    cout << "Results saved to " << cl_params["output"] << endl;
    for(int i = 1; i < argc; i++){
      if (string(argv[i]).compare("-d") == 0 || string(argv[i]).compare("-dump") == 0){
        dump_matrix(matrix, cl_params["dump"]);
        cout << "Matrix saved to " << cl_params["dump"] << endl;
      }
    }
    return 0;
}