#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <unordered_map>
#include <bitset>
#include <queue>
#include <filesystem>
using namespace std;
const int char_buffered = 2147483647;
const int window_size = 16000;
const int Buffer_size = 128;
//1048576

class Token {
public:
    int size;
    int length;
    char next;
    bool flag_next;
};

class LZ78_codes {
public:
    int index;
    short int next_сhar;
};

class element_tree {
public:
    char symbol;
    int count;
    element_tree* left_child = nullptr;
    element_tree* right_child = nullptr;
    element_tree(char new_symbol, int new_count){
        symbol = new_symbol;
        count = new_count;
        left_child = nullptr;
        right_child = nullptr;
    }
};

class compare {
public:
    bool operator()(element_tree* left, element_tree* right) {
        return left->count > right->count;
    }
};

void read_file_by_bytes(const string& filePath, vector <short int>& mass) 
{
    ifstream file(filePath, ios::binary);
    if (!file.is_open()) 
    {
        cout << "Не удалось открыть файл: " << filePath << endl;
        return;
    }

    char byte;
    while (file.get(byte)) 
    {
        mass.push_back(static_cast<short int>(byte));
    }
    file.close();
}

void write_by_bytes(const string& filePatch, vector<short int>& mass) {
    ofstream file(filePatch, ios::binary);
    if (!file.is_open())
    {
        cout << "Не удалось открыть файл: " << filePatch << endl;
        return;
    }
    for (int i = 0; i < mass.size(); i++) {
        file.write(reinterpret_cast<const char*>(&mass[i]), 1);
    }
    file.close();
}

void shell_sort(vector<vector<short int>>& matrix) {
    int n = matrix.size();
    int gap = n / 2;
    while (gap > 0) {
        for (int i = 0; i < gap; i++) {
            int j = gap;
            while (j < n) {
                while (j >= gap and matrix[j] < matrix[j - gap]) {
                    /*matrix = swap(matrix, j, j - gap);*/
                    vector<short int> buffer = matrix.at(j);
                    matrix.at(j) = matrix.at(j - gap);
                    matrix.at(j - gap) = buffer;
                    j -= gap;
                }
                j += gap;
            }
        }
        gap /= 2;
    }
}

void rle_compress(vector<short int> str, vector<short int>& new_str) {
    int index = 0;
    int len = str.size();
    bool reapet;
    vector<short int> buffer_not_reapet;
    vector<short int> buffer_reapet;
    if (str[0] == str[1]) {
        new_str.push_back(0);
        buffer_reapet.push_back(str[index]);
        index++;
        buffer_reapet.push_back(str[index]);
        index++;
        reapet = true;
    }
    else {
        buffer_not_reapet.push_back(str[index]);
        index++;
        buffer_not_reapet.push_back(str[index]);
        index++;
        reapet = false;
    }
    while (index < len - 1) {
        if (reapet) {
            if (str[index - 1] == str[index]) {
                buffer_reapet.push_back(str[index]);
                if (buffer_reapet.size() == 127) {
                    new_str.push_back(buffer_reapet.size());
                    new_str.push_back(buffer_reapet[0]);
                    buffer_reapet.clear();
                    new_str.push_back(0);
                }
            }
            else {
                buffer_not_reapet.push_back(str[index]);
                reapet = false;
                if (buffer_reapet.size() > 0) {
                    new_str.push_back(buffer_reapet.size());
                    new_str.push_back(buffer_reapet[0]);
                }
                buffer_reapet.clear();
                if (new_str[new_str.size() - 1] == 0 and str[index - 1] != 0) {
                    new_str.pop_back();
                }
            }
        }
        else {
            if (str[index - 1] != str[index]) {
                buffer_not_reapet.push_back(str[index]);
                if (buffer_not_reapet.size() == 127) {
                    if (str[index + 1] == str[index]) {
                        buffer_not_reapet.pop_back();
                    }
                    new_str.push_back(buffer_not_reapet.size());
                    for (int i = 0; i < buffer_not_reapet.size(); i++) {
                        new_str.push_back(buffer_not_reapet[i]);
                    }
                    buffer_not_reapet.clear();
                    new_str.push_back(0);
                    new_str.push_back(0);
                }
            }
            else {
                if (buffer_not_reapet.size() > 0) {
                    buffer_not_reapet.pop_back();
                }
                buffer_reapet.push_back(str[index - 1]);
                buffer_reapet.push_back(str[index]);
                new_str.push_back(buffer_not_reapet.size());
                for (int i = 0; i < buffer_not_reapet.size(); i++) {
                    new_str.push_back(buffer_not_reapet[i]);
                }
                buffer_not_reapet.clear();
                reapet = true;
            }
        }
        index++;
    }
    if (reapet) {
        buffer_reapet.push_back(str[index]);
        new_str.push_back(buffer_reapet.size());
        new_str.push_back(buffer_reapet[0]);
    }
    else {
        buffer_not_reapet.push_back(str[index]);
        new_str.push_back(buffer_not_reapet.size());
        for (int i = 0; i < buffer_not_reapet.size(); i++) {
            new_str.push_back(buffer_not_reapet[i]);
        }
    }
}

void rle_decompress(vector<short int> str, vector<short int>& old_str) {
    int index = 0;
    int len = str.size();
    bool reapet = false;
    int len_input_non_reapet = str[index];
    int len_input_reapet = 0;
    short int copy_element = 0;
    index++;
    while (index < len) {
        if (not reapet) {
            for (int i = 0; i < len_input_non_reapet; i++) {
                old_str.push_back(str[index]);
                index++;
            }
            reapet = true;
            if (index < len) {
                len_input_reapet = str[index];
            }
            index++;
            if (index < len) {
                copy_element = str[index];
            }
        }
        else {
            for (int i = 0; i < len_input_reapet; i++) {
                old_str.push_back(str[index]);
            }
            index++;
            reapet = false;
            if (index < len) {
                len_input_non_reapet = str[index];
            }
            index++;
        }
    }
}

void BWT_code(vector<short int> str, vector<short int>& bwt_str, short int& position_bwt) {
    vector<vector<short int>> bwt_mass;
    bwt_mass.push_back(str);
    for (int i = 1; i < str.size(); i++) {
        vector<short int> new_str;
        for (int j = 1; j < bwt_mass[i - 1].size(); j++) {
            new_str.push_back(bwt_mass[i - 1][j]);
        }
        new_str.push_back(bwt_mass[i - 1][0]);
        bwt_mass.push_back(new_str);
    }
    shell_sort(bwt_mass);
    for (int i = 0; i < bwt_mass.size(); i++) {
        if (bwt_mass[i] == str) {
            position_bwt = i;
        }
        bwt_str.push_back(bwt_mass[i][bwt_mass[i].size() - 1]);
    }
}

void buffer_for_BWT_code(vector<short int> str, vector<short int>& bwt_str, vector<short int>& position_bwt) {
    int len_buffer = 1200;
    int n = str.size();
    int len;
    long long index = 0;
    vector<short int> buffer;
    vector<short int> bwt_buffer;
    short int position_bwt_ = 0;
    while (index < n) {
        buffer.push_back(str[index]);
        if (buffer.size() % len_buffer == 0) {
            BWT_code(buffer, bwt_buffer, position_bwt_);
            len = bwt_buffer.size();
            for (int i = 0; i < len; i++) {
                bwt_str.push_back(bwt_buffer[i]);
            }
            position_bwt.push_back(position_bwt_);
            position_bwt_ = 0;
            buffer.clear();
            bwt_buffer.clear();
        }
        index++;
    }
    if (buffer.size() > 0) {
        BWT_code(buffer, bwt_buffer, position_bwt_);
        len = bwt_buffer.size();
        for (int i = 0; i < len; i++) {
            bwt_str.push_back(bwt_buffer[i]);
        }
        position_bwt.push_back(position_bwt_);
        position_bwt_ = 0;
        buffer.clear();
        bwt_buffer.clear();
    }
}

void count_sort(vector<short int>& P_inverse, vector<short int> str) {
    int n = 256;
    int len = str.size();
    vector<short int> T;
    vector<short int> T_sub;
    for (int i = 0; i < n; i++) {
        T.push_back(0);
        T_sub.push_back(0);
    }
    for (int i = 0; i < len; i++) {
        T[str[i] + 128]++;
    }
    for (int i = 1; i < n; i++) {
        T_sub[i] = T_sub[i - 1] + T[i - 1];
    }
    for (int i = 0; i < len; i++) {
        P_inverse.push_back(-1);
    }
    for (int i = 0; i < len; i++) {
        P_inverse[T_sub[str[i] + 128]] = i;
        T_sub[str[i] + 128]++;
    }
}

void BWT_decode(vector<short int> bwt_str, vector<short int>& str, int position_bwt) {
    int len = bwt_str.size();
    vector<short int>P_inverse;
    count_sort(P_inverse, bwt_str);

    int j = position_bwt;
    for (int i = 0; i < len; i++) {
        j = P_inverse[j];
        str.push_back(bwt_str[j]);
    }
}

void buffer_for_BWT_decode(vector<short int> bwt_str, vector<short int>& str, vector<short int>& position_bwt) {
    int len_buffer = 1200;
    int n = bwt_str.size();
    int len;
    long long index = 0;
    long long bufers_index = 0;
    vector<short int> buffer;
    vector<short int> bwt_buffer;
    while (index < n) {
        buffer.push_back(bwt_str[index]);
        if (buffer.size() == len_buffer) {
            BWT_decode(buffer, bwt_buffer, position_bwt[bufers_index]);
            len = bwt_buffer.size();
            for (int i = 0; i < len; i++) {
                str.push_back(bwt_buffer[i]);
            }
            buffer.clear();
            bwt_buffer.clear();
            bufers_index++;
        }
        index++;
    }
    if (buffer.size() > 0) {
        BWT_decode(buffer, bwt_buffer, position_bwt[bufers_index]);
        len = bwt_buffer.size();
        for (int i = 0; i < len; i++) {
            str.push_back(bwt_buffer[i]);
        }
        buffer.clear();
        bwt_buffer.clear();
        bufers_index++;
    }
}

void count_symbols_for_HA(string input, unordered_map<char, int>& counter) {
    ifstream input_file(input, ios::binary);
    const size_t buffer_size = char_buffered;
    vector<char> buffer(buffer_size);
    while (input_file) {
        input_file.read(buffer.data(), buffer_size);
        size_t bytes_readed = input_file.gcount();
        for (size_t i = 0; i < bytes_readed; i++) {
            counter[buffer[i]]++;
        }
    }
}

element_tree* HA_building_tree(unordered_map<char, int>& counter) {
    priority_queue<element_tree*, vector<element_tree*>, compare> Heap;

    for (auto pair : counter) {
        Heap.push(new element_tree(pair.first, pair.second));
    }
    while (Heap.size() > 1) {
        element_tree* left = Heap.top(); Heap.pop();
        element_tree* right = Heap.top(); Heap.pop();
        element_tree* new_element = new element_tree('$', left->count + right->count);
        new_element->left_child = left;
        new_element->right_child = right;
        Heap.push(new_element);
    }
    return (Heap.top());
}

void write_codes(element_tree* root, string str, unordered_map<char, string>& huffmanCodes) {
    if (!root) {
        return;
    }

    if (root->left_child == nullptr and root->right_child == nullptr) {
        cout << root->symbol << "-" << root->count << ":" << str << endl;
        huffmanCodes[root->symbol] = str;
    }

    write_codes(root->left_child, str + "0", huffmanCodes);
    write_codes(root->right_child, str + "1", huffmanCodes);
}

void compressed_HA(const string& input, const string& output) {
    unordered_map<char, int> counter;
    unordered_map<char, string> huffman_codes;
    element_tree* root = nullptr;
    ifstream input_file(input, ios::binary);
    ofstream file(output, ios::binary);

    count_symbols_for_HA(input, counter);
    root = HA_building_tree(counter);
    write_codes(root, "", huffman_codes);
    size_t codes_size = counter.size();
    cout << "s" << endl;
    file.write(reinterpret_cast<const char*>(&codes_size), sizeof(codes_size));
    for (const auto& part : counter) {
        file.write(&part.first, sizeof(part.first));
        file.write(reinterpret_cast<const char*>(&part.second), sizeof(part.second));
    }
    cout << "s" << endl;
    vector<char> input_buffer(sizeof(input_file));
    string string_of_bits;

    while (input_file) {
        input_file.read(input_buffer.data(), sizeof(input_file));
        size_t bytes_readed = input_file.gcount();

        for (size_t i = 0; i < bytes_readed; i++) {
            string_of_bits += huffman_codes[input_buffer[i]];
        }
        cout << "s";
        while (string_of_bits.size() >= 8) {
            bitset<8> bits(string_of_bits.substr(0, 8));
            char byte = static_cast<char>(bits.to_ulong());
            file.write(&byte, 1);
            string_of_bits.erase(0, 8);
        }
    }
    cout << "s" << endl;
    if (!string_of_bits.empty()) {
        while (string_of_bits.size() < 8) {
            string_of_bits = string_of_bits + '0';
        }
        bitset<8> bits(string_of_bits);
        char byte = static_cast<char>(bits.to_ulong());
        file.write(&byte, 1);
    }

    input_file.close();
    file.close();
}

void decompressed_HA(const string& input, const string& output) {
    ifstream compressed(input, ios::binary);
    ofstream decompressed(output, ios::binary);
    size_t tableSize;
    compressed.read(reinterpret_cast<char*>(&tableSize), sizeof(tableSize));
    unordered_map<char, int> counter;
    element_tree* root = nullptr;
    unsigned long long count_symbols = 0;

    for (size_t i = 0; i < tableSize; i++) {
        char buffer_char;
        int count;
        compressed.read(&buffer_char, sizeof(buffer_char));
        compressed.read(reinterpret_cast<char*>(&count), sizeof(count));
        cout << endl << count;

        counter[buffer_char] = count;
        cout << count << endl;
    }
    for (const auto& part : counter) {
        count_symbols += part.second;
    }
    root = HA_building_tree(counter);

    vector<char> buffer(sizeof(compressed));
    string string_of_bits;
    element_tree* this_child = root;
    while (compressed) {
        compressed.read(buffer.data(), sizeof(compressed));
        size_t bytes_readed = compressed.gcount();

        for (size_t i = 0; i < bytes_readed; i++) {
            bitset<8> bits(buffer[i]);
            string_of_bits += bits.to_string();
        }

        size_t pos = 0;
        while (pos < string_of_bits.size()) {
            if (string_of_bits[pos] == '0') {
                this_child = this_child->left_child;
            }
            else {
                this_child = this_child->right_child;
            }

            if ((this_child->left_child == nullptr) and (this_child->right_child == nullptr)) {
                if (count_symbols > 0) {
                    decompressed.write(&this_child->symbol, 1);
                    count_symbols--;
                }
                this_child = root;
            }
            pos++;
        }

        if (pos > 0) {
            string_of_bits.erase(0, pos);
        }
    }

    compressed.close();
    decompressed.close();
}

void MTF_code(vector<short int> str, vector<short int>& new_str) {
    int len = str.size();
    int index;
    vector<short int> shifts;
    for (int i = -128; i < 128; i++) {
        shifts.push_back(static_cast<short int>(static_cast<char>(i)));
    }
    //cout << "s" << endl;
    for (int i = 0; i < len; i++) {
        index = 0;
        //cout << shifts.size() << "|" << str.size() << endl;
        //cout << str[i] << "|" << shifts[index] << endl;
        while (shifts[index] != str[i]) {
            //cout << str[i] << "|" << shifts[index] << endl;
            index++;
        }
        //cout << "s";
        new_str.push_back(index - 128);
        vector<short int> buffer;
        buffer.push_back(shifts[index]);
        for (int j = 0; j < index; j++) {
            buffer.push_back(shifts[j]);
        }
        for (int j = index + 1; j < shifts.size(); j++) {
            buffer.push_back(shifts[j]);
        }
        shifts = buffer;
        //cout << i << " :" << str[i] << endl;
    }
}

void MTF_decode(vector<short int> str, vector<short int>& old_str) {
    int len = str.size();
    int index;
    vector<short int> shifts;
    for (int i = -128; i < 128; i++) {
        shifts.push_back(static_cast<short int>(static_cast<char>(i)));
    }
    for (int i = 0; i < len; i++) {
        index = static_cast<char>(str[i]) + 128;
        old_str.push_back(shifts[index]);
        vector<short int> buffer;
        buffer.push_back(shifts[index]);
        for (int j = 0; j < index; j++) {
            buffer.push_back(shifts[j]);
        }
        for (int j = index + 1; j < shifts.size(); j++) {
            buffer.push_back(shifts[j]);
        }
        shifts = buffer;

    }
}

void lz78_compress(const string& input, const string& output) {
    ifstream input_file(input, ios::binary);
    vector<unsigned char> readed_data((istreambuf_iterator<char>(input_file)), istreambuf_iterator<char>());
    input_file.close();
    vector<LZ78_codes> encoded_readed_data;
    unordered_map<string, int> dictionary;
    string buffer;
    string write_buffer;
    int next_index = 1;
    int index;

    for (unsigned char symbol : readed_data) {
        buffer.push_back(symbol);
        if (dictionary.find(buffer) == dictionary.end()) {
            dictionary[buffer] = next_index;
            next_index++;
            write_buffer = buffer;
            write_buffer.pop_back();
            index = (buffer.size() == 1) ? 0 : dictionary[write_buffer];
            encoded_readed_data.push_back({ index, symbol });
            buffer.clear();
        }
    }

    if (!buffer.empty()) {
        write_buffer = buffer;
        write_buffer.pop_back();
        index = dictionary[write_buffer];
        encoded_readed_data.push_back({ index, buffer.back() });
    }

    ofstream output_file(output, ios::binary);

    for (const auto& encoded : encoded_readed_data) {
        output_file.write(reinterpret_cast<const char*>(&encoded.index), sizeof(encoded.index));
        output_file.write(reinterpret_cast<const char*>(&encoded.next_сhar), sizeof(encoded.next_сhar));
    }

    output_file.close();
}

void lz78_decompress(const string& input, const string& output) {
    ifstream input_file(input, ios::binary);
    vector<LZ78_codes> encoded_data;
    LZ78_codes encoded;

    while (input_file.read(reinterpret_cast<char*>(&encoded.index), sizeof(encoded.index)) and input_file.read(reinterpret_cast<char*>(&encoded.next_сhar), sizeof(encoded.next_сhar))) {
        encoded_data.push_back(encoded);
    }
    input_file.close();

    vector<unsigned char> decoded_data;
    unordered_map<int, string> dictionary;
    int next_index = 1;

    for (const auto& encoded : encoded_data) {
        int index = encoded.index;
        unsigned char next_сhar = encoded.next_сhar;

        string entry = dictionary[index];
        entry.push_back(next_сhar);
        dictionary[next_index++] = entry;

        decoded_data.insert(decoded_data.end(), entry.begin(), entry.end());
    }

    ofstream output_file(output, ios::binary);
    output_file.write(reinterpret_cast<const char*>(decoded_data.data()), decoded_data.size());
    output_file.close();
}

void lz77_compress(const string& input_file, const string& output) {
    vector <char> matrix;
    char byte;
    ifstream file(input_file, ios::binary);
    while (file.get(byte)) {
        matrix.push_back(byte);
    }
    file.close();
    vector<Token> tokens;
    long long i = 0;
    int len_matrix = matrix.size();
    while (i < len_matrix) {
        int bestlen = 0;
        int size = 0;
        int index = 0;
        if (i >= window_size) {
            index -= window_size;
        }
        for (int j = index; j < i; j++) {
            int length = 0;
            while (length < Buffer_size and i + length < matrix.size() and matrix[j + length] == matrix[i + length]) {
                length++;
            }
            if (length > bestlen) {
                bestlen = length;
                size = i - j;
            }
        }
        Token token;
        if (i + bestlen < matrix.size()) {

            token.size = (bestlen > 0) ? size : 0;
            token.length = bestlen;
            token.next = matrix[i + bestlen];
            token.flag_next = false;
            tokens.push_back(token);
            i += bestlen + 1;
        }
        else {
            token.size = (bestlen > 0) ? size : 0;
            token.length = bestlen;
            token.next = 0;
            token.flag_next = true;
            tokens.push_back(token);
            i += bestlen;
        }
    }


    ofstream output_file(output, ios::binary);
    unsigned int tokenCount = tokens.size();

    output_file.write(reinterpret_cast<const char*>(&tokenCount), sizeof(tokenCount));
    for (const auto& token : tokens) {
        output_file.write(reinterpret_cast<const char*>(&token.size), sizeof(token.size));
        output_file.write(reinterpret_cast<const char*>(&token.length), sizeof(token.length));
        output_file.write(&token.next, sizeof(token.next));
        output_file.write(reinterpret_cast<const char*>(&token.flag_next), sizeof(token.flag_next));
    }

    output_file.close();
}

void lz77_decompress(const string& input, const string& output) {
    ifstream input_file(input, ios::binary);
    int tokenCount;
    input_file.read(reinterpret_cast<char*>(&tokenCount), sizeof(tokenCount));
    vector<Token> tokens(tokenCount);

    for (int i = 0; i < tokenCount; i++) {
        Token token;
        input_file.read(reinterpret_cast<char*>(&token.size), sizeof(token.size));
        input_file.read(reinterpret_cast<char*>(&token.length), sizeof(token.length));
        input_file.read(&token.next, sizeof(token.next));
        input_file.read(reinterpret_cast<char*>(&token.flag_next), sizeof(token.flag_next));
        tokens[i] = token;
    }
    input_file.close();

    vector <char> matrix;
    for (const auto& token : tokens) {
        if (token.length > 0) {
            int pos = matrix.size() - token.size;
            for (int i = 0; i < token.length; i++) {
                matrix.push_back(matrix[pos + i]);
            }
        }
        if (not token.flag_next) {
            matrix.push_back(token.next);
        }
    }

    ofstream output_file(output, ios::binary);

    output_file.write(matrix.data(), matrix.size());
    output_file.close();
}

/*
DTLiteInstaller-11.exe 
enwik7.txt
martyr.txt
*/
/*
exe_compress.txt
enwik_compress.txt
ru_text_compress.txt
*/
/*
BW.raw
BWG.raw
RGB.raw

test.txt
untest.txt
*/
int main()
{
    vector <short int> str;
    vector <short int> str_1;
    vector <short int> empty;
    float szhat;
    string str_f;
    szhat;
    /*for (int windowSize = 128; windowSize <= 32768; windowSize *= 2) {
        lz77_compress("TEST_FILES\\enwik7.txt", "TEST_GRAPHICS\\enwik_compress.txt", windowSize);
    }*/
    //lz78_compress("TEST_FILES\\enwik7.txt", "test.txt");
    //lz78_compress("test.txt", "untest.txt");
    //lz78_decompress("test.txt", "untest.txt");
    //lz78_decompress("LZ78\\enwik_compress.txt", "untest.txt");
}