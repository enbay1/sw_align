#include <iostream>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <filesystem>
using namespace std;

vector<vector<string>> read_files(vector<string> files){
    vector<vector<string>> return_vector;
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
                file_contents.push_back(line);
            }
            input_file.close();
        }
        return_vector.push_back(file_contents);
    }
    return return_vector;
}

map<string, string> parse_cl(int argc, char **argv){
    map<string, string> parameters;
    // Program needs at least <program name> <file1> <file2> in argv
    if (argc < 3){
        cout << "At least two input files need to be specified" << endl;
    }
    // Check that the files are of type .fasta
    for (string file: {argv[1], argv[2]}){
        if (!file.find(".fasta") && !file.find(".FASTA")){
            cerr << "File "  + file + " not of type .fasta";
            exit(1);
        }
    }
    parameters.emplace("file_1", argv[1]);
    parameters.emplace("file_2", argv[2]);
    for(int i = 3; i < argc; i++){
        cout << "Flag found" << endl;
    }
    return parameters;
}

int main(int argc, char **argv){
    map<string, string> cl_params;
    vector<vector<string>> file_contents;
    cl_params = parse_cl(argc, argv);
    vector<string> passable_file_names;
    passable_file_names.push_back(cl_params["file_1"]);
    passable_file_names.push_back(cl_params["file_2"]);
    file_contents = read_files(passable_file_names);
    return 0;
}

/*    
    vector<string> v;
    v.push_back("one");
    for ( auto &i : v ) {
        cout << i << endl;
    }
    */