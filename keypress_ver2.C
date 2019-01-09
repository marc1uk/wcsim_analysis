/* fakeKey.c */
// from https://bharathisubramanian.wordpress.com/2010/03/14/x11-fake-key-event-generation-using-xtest-ext/
#include <X11/Xlib.h>
#include <X11/Intrinsic.h>
#include <X11/extensions/XTest.h>
#include <unistd.h>
/* Send Fake Key Event */
static void SendKey (Display * disp, KeySym keysym, KeySym modsym){
 KeyCode keycode = 0, modcode = 0;
 keycode = XKeysymToKeycode (disp, keysym);
 if (keycode == 0) return;
 XTestGrabControl (disp, True);
 /* Generate modkey press */
 if (modsym != 0) {
  modcode = XKeysymToKeycode(disp, modsym);
  XTestFakeKeyEvent (disp, modcode, True, 0);
 }
 /* Generate regular key press and release */
 XTestFakeKeyEvent (disp, keycode, True, 0);
 XTestFakeKeyEvent (disp, keycode, False, 0); 
 
 /* Generate modkey release */
 if (modsym != 0)
  XTestFakeKeyEvent (disp, modcode, False, 0);
 
 XSync (disp, False);
 XTestGrabControl (disp, False);
}
 
/* Main Function */
int main (){
 Display *disp = XOpenDisplay (NULL);
 sleep (5);
 /* Send ASCII A & B */
 SendKey (disp, XK_A, 0);
 SendKey (disp, XK_B, 0);
 /* Send ALT+Tab */
 sleep (3);
 SendKey (disp, XK_Tab, XK_Alt_L);
 sleep (3);
 SendKey (disp, XK_Tab, XK_Alt_L);
}

Build and Run it using following commands:
$ gcc fakeKey.c -o fakeKey -lX11 -lXtst -lXext
$ ./fakeKey

Look into /usr/include/X11/ keysymdef.h header file to know about key symbols of other keys.
