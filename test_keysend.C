// Send a fake keystroke event to an X window.
// by Adam Pierce - http://www.doctort.org/adam/
// This is public domain software. It is free to use by anyone for any purpose.
// compile with: `g++ -std=c++11 test_keysend.C -o testkeysend -L/usr/lib/x86_64-linux-gnu -lX11 -lpthread -lXtst`
// requires the libxtst-dev package to be installed

#include <X11/Xlib.h>
#include <X11/keysym.h>
#include <X11/extensions/XTest.h> // /usr/include/X11/extensions/XTest.h

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <thread>   // thread
#include <future>   // promise
//#include <unistd.h> // usleep

// The key code to be sent.
// A full list of available codes can be found in /usr/include/X11/keysymdef.h
//#define KEYCODE XK_Down

void ExternalFunction(std::promise<int> finishedin);

// Function to create a keyboard event
XKeyEvent createKeyEvent(Display *display, Window &win,
                           Window &winRoot, bool press,
                           int keycode, int modifiers)
{
   XKeyEvent event;
   
   event.display     = display;
   event.window      = win;
   event.root        = winRoot;
   event.subwindow   = None;
   event.time        = CurrentTime;
   event.x           = 1;
   event.y           = 1;
   event.x_root      = 1;
   event.y_root      = 1;
   event.same_screen = True;
   event.keycode     = XKeysymToKeycode(display, keycode);
   event.state       = modifiers;
   
   if(press)
      event.type = KeyPress;
   else
      event.type = KeyRelease;
   
   return event;
}


int main(int argc, char* argv[]){
	
	// Obtain the X11 display.
	Display *display = XOpenDisplay(0);
	if(display == NULL)
	return -1;
	// Note multiple calls to XOpenDisplay() with the same parameter return different handles;
	// sending the press event with one handle and the release event with another will not work.
	
	// Get the root window for the current display.
	Window winRoot = XDefaultRootWindow(display);
	
	// We need to prevent auto-repeat in order to send single keypresses, 
	// otherwise the keydown and keyup commands result in multiple presses
	XAutoRepeatOff(display);
	//XChangeKeyboardControl   // for google reference
	//XAutoRepeatOn(display);  // for reference, if needed
	// XXX XXX XXX THIS AFFECTS THE GLOBAL ENVIRONMENT - YOU NEED TO RE-ENABLE IT
	// OR YOUR COMPUTER WILL STOP REPEATING KEYS AFTER THE PROGRAM CLOSES!!!! XXX 
	
	
	// Start the program we want to receieve the keypress
	// Start in an external thread so it can execute while we continue and send keys to it
	// first create a promise we can use to hold until the external thread is done
	std::promise<int> barrier;
	std::future<int> barrier_future = barrier.get_future();
	
	// start the thread
	std::cout<<"starting thread"<<std::endl;
	std::thread athread(ExternalFunction, std::move(barrier));
	
	// wait for just a little bit for it to start
	std::this_thread::sleep_for(std::chrono::seconds(2));
	
	// Find the window which has the current keyboard focus.
	Window winFocus;
	int    revert;
	XGetInputFocus(display, &winFocus, &revert);
	
	// Send a fake key press event to the window.
	std::cout<<"sending keydown"<<std::endl;
	int KEYCODE = XK_S;  // XK_s sends lower case s etc.
	XKeyEvent event = createKeyEvent(display, winFocus, winRoot, true, KEYCODE, 0);
	//XSendEvent(event.display, event.window, True, KeyPressMask, (XEvent *)&event);
	// This may not work for sending key events to the focused window, because
	// key events generated in this way are tagged as fake and many programs 
	// ignore events with that tag. So use the XTest function instead.
	// www.handhelds.org/moin/moin.cgi/GeneratingSyntheticX11Events
	XTestFakeKeyEvent(event.display, event.keycode, True, CurrentTime);
	
	// You may need to request the focus window for every new key event you send, i.e.
	// getfocus(), keypress(), keyrelease() will not work
	// getfocus(), keypress(), getfocus(), keyrelease() will work
	XGetInputFocus(display, &winFocus, &revert);
	
	//std::cout<<"sending key release"<<std::endl;
	// Send a fake key release event to the window.
	//event = createKeyEvent(display, winFocus, winRoot, false, KEYCODE, 0); // << doesn't work/not needed
	//XSendEvent(event.display, event.window, True, KeyPressMask, (XEvent *)&event);
	XTestFakeKeyEvent(event.display, event.keycode, False, CurrentTime);
	// NOTE : FALSE to release the key: changing createKeyEvent doesn't seem to work
	
	// wait for the external thread to indicate completion
	std::cout<<"waiting for external thread to complete"<<std::endl;
	std::chrono::milliseconds span(100);
	std::future_status thestatus;
	do {
		thestatus = barrier_future.wait_for(span);
	} while (thestatus==std::future_status::timeout);
	// get the return to know if the external thead succeeded
	std::cout << "thread returned " << barrier_future.get() << std::endl;
	
	// call join to cleanup the external thread
	athread.join();
	
	// XXX RE-ENABLE AUTO-REPEAT! XXX
	XAutoRepeatOn(display);
	
	// Cleanup X11 display
	XCloseDisplay(display);
	
	return 1;
}

void ExternalFunction(std::promise<int> finishedin){
	
	// if we need to pass between member functions in classes, must use std::move
	// otherwise we can just call set_value on finishedin
	std::promise<int> finished;
	finished = std::promise<int>(std::move(finishedin));
	
	// start the external program, which will wait for a keypress then return
	const char* command = "./kbhit_test";
	std::cout<<"calling command "<<command<<std::endl;
	system(command);
	std::cout<<"returned from command"<<std::endl;
	
	// release the barrier, allowing the caller to continue
	finished.set_value(1);
}
