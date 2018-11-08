#include <iostream>
#include <thread>         // std::thread
#include <mutex>          // std::mutex, std::lock_guard
#include <stdexcept>      // std::logic_error

//std::mutex mtx;

class myclass{
  public:
  myclass(){};
  ~myclass(){};
  int something;
  void DoStuff(std::mutex* mtx){
    std::lock_guard<std::mutex> lck (*mtx); // lock the mutex. Will unlock when the function ends.
    for(int i=0; i<3; i++){
      std::cout << i << " ... ";
      std::this_thread::sleep_for (std::chrono::seconds(1));
    }
    std::cout<<std::endl;
  }
};

int main(){
//  std::thread threads[10];
//  // spawn 10 threads:
//  for (int i=0; i<10; ++i)
//    threads[i] = std::thread(print_thread_id,i+1);
//  for (auto& th : threads) th.join();
  
  std::mutex* mtx = new std::mutex();
  myclass aclass;
  aclass.DoStuff(mtx);
  std::lock_guard<std::mutex> lck (*mtx); // lock the mutex. Requires the mutex to be free.
  std::cout<<"function has ended"<<std::endl;
  
  return 0;
}
