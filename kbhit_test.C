// C++ program to demonstrate use of kbhit() 
#include <iostream>
//#include <conio>        // kbhit, for c. windows only
#include "keyb.h"         // linux ver, from CAENs
#include <thread>         // std::this_thread::sleep_for
#include <chrono>         // std::chrono::seconds
//#include <unistd.h>     // usleep, if needed to replace sleep_for
// compile with g++ -std=c++11 kbhit_test.C keyb.c -o kbhit_test

int main()
{ 
    while (!kbhit()){
        printf("Press a key\n");
        std::this_thread::sleep_for (std::chrono::seconds(1));
        //usleep(1000000);  //microseconds
    }
    return 0;
}
