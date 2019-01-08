#include <iostream>
#include <thread>         // std::thread
#include <future>
#include <stdexcept>      // std::logic_error
// must compile with  -pthread option

class myclass{
  public:
  myclass(){};
  ~myclass(){};
  int something;
  void DoStuff(std::promise<void> barrier){
    for(int i=0; i<3; i++){
      std::cout << i << " ... ";
      std::this_thread::sleep_for (std::chrono::seconds(1));
    }
    std::cout<<std::endl;
    barrier.set_value();
  }
};

int main(){
//  std::thread threads[10];
//  // spawn 10 threads:
//  for (int i=0; i<10; ++i)
//    threads[i] = std::thread(print_thread_id,i+1);
//  for (auto& th : threads) th.join();
  
  std::promise<void> barrier;  // will hold until we're done viewing waveforms
  std::future<void> barrier_future = barrier.get_future();
  myclass aclass;
  aclass.DoStuff(std::move(barrier));
  barrier_future.wait();
  std::cout<<"function has ended"<<std::endl;
  
  return 0;
}
