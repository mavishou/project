// #include <string>
// #include <iostream>
// #include <stdio.h>
// using namespace std;

// #include<stdio.h> 
        
// int main() 
        
// { 
        
//     FILE * fp; 
//     char buffer[80]; 
//     fp=popen("readlink -f $0","r"); 
//     fgets(buffer, sizeof(buffer),fp); 
//     printf("%s",buffer); 
//     pclose(fp); 
//     return 0;    
// } 


#include <iostream>
#include <ostream>
#include <boost/filesystem/convenience.hpp>

int main(int argc, char** argv)
{
    boost::filesystem::path argvPath( argv[0] );
    boost::filesystem::path executablePath( argvPath.parent_path() );
    boost::filesystem::path runPath( boost::filesystem::initial_path() );

    std::cout << executablePath << std::endl;
    std::cout << runPath << std::endl;
    return 0;
}