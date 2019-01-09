//#include <string>
#include <iostream>
//#include <sstream>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char* argv[]){
	
	const char* command = "./kbhit_test";
	std::cout<<"calling command "<<command<<std::endl;
	system(command);
	std::cout<<"returned from command"<<std::endl;
	
//	// use this to read the stdout return.
//	FILE* file = popen(command, "r");
//	
//	std::string command_output="";
//	if(file){
//		std::ostringstream stm ;
//		constexpr std::size_t MAX_LINE_SZ = 1024 ;
//		char line[MAX_LINE_SZ] ;
//		while( fgets( line, MAX_LINE_SZ, file ) ) stm << line << '\n' ;
//		pclose(file) ;
//		command_output = stm.str();
//	}
//	
	return 1;
}
