// This is a dumb code where I test random stuff that doesn't work in other codes

#include <iostream>
#include <string>

int main(int argc, char** argv)
{
    std::string str = "/path/to/file.root";

    std::cout << "Splitting: " << str << '\n';
    std::size_t found = str.find_last_of("/\\");
    std::cout << " path: " << str.substr(0,found) << '\n';
    std::cout << " file: " << str.substr(found+1) << '\n';

}
