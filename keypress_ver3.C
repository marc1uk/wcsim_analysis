/* genxkey - small example program to generate X11 keyboard events

   (c) 2005 Robin Gareus <robin@gareus.org>
   
   Compile with: 

     gcc -o genxkey xkbd.c -lX11 -lXtst 
  
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software Foundation,
   Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.  

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include <X11/keysym.h>
#include <X11/keysymdef.h>

#include <X11/extensions/XTest.h>

Display   *display = NULL;

void usage(int exit_status, char *program_name) {
	fprintf(stderr,"Send a keyboard event to the currently focused X-application.\n",program_name);
	fprintf(stderr,"(c) GPL 2005 Robin Gareus <robin@gareus.org>\n");
	fprintf(stderr," usage: %s [-d <int>] [-p] <key>\n",program_name);
	fprintf(stderr,"   -p : do not actually send the event, but print the KeySym info.\n");
	fprintf(stderr,"   -d <int>: delay event by <int> seconds\n");
	fprintf(stderr," Standard KeySym names are obtained from <X11/keysymdef.h> by removing the\n XK_ prefix from each name. Each KeySym may be modified by prepending \n either of ALT-, META-, CTRL- or SHIFT-.\n");
	exit(exit_status);
}

int xkey_init (char *display_name) {

	if (display_name == NULL) display_name = XDisplayName (NULL);
//	printf ("Display name: %s\n", display_name);

  	if ( (display=XOpenDisplay(display_name)) == NULL ) {
		fprintf( stdout, "Cannot connect to X server\n");
		return(-1);
	}
	return (0);
}

	//ks = XStringToKeysym("w");
	//mod = 0: no modifier - 1: Meta/Alt key pressed. - 2: CTRL key pressed
void xkey_send (long ks, int mod) {
	int kc;

	if (mod) { 
		if (mod==1) kc = XKeysymToKeycode (display, XK_Meta_L);
		else if (mod==3) kc = XKeysymToKeycode (display, XK_Alt_L);
		else if (mod==4) kc = XKeysymToKeycode (display, XK_Shift_L);
		else kc = XKeysymToKeycode (display, XK_Control_L);
		XTestFakeKeyEvent(display, kc, True, 0);
		XFlush(display);
	}

	kc = XKeysymToKeycode (display, ks);
//	fprintf(stdout, "XTestFakeKeyEvent(0x%lx, %ld, %d)\n", (long)display, (long)kc, type == KeyPress);


	XTestFakeKeyEvent(display, kc, True, 0);
	XFlush(display);
	XTestFakeKeyEvent(display, kc, False, 0);
	XFlush(display);

	if (mod) { 
		if (mod==1) kc = XKeysymToKeycode (display, XK_Meta_L);
		else if (mod==3) kc = XKeysymToKeycode (display, XK_Alt_L);
		else if (mod==4) kc = XKeysymToKeycode (display, XK_Shift_L);
		else kc = XKeysymToKeycode (display, XK_Control_L);
		XTestFakeKeyEvent(display, kc, False, 0);
		XFlush(display);
	}
}

const char *modifier[] = { "none", "META", "CTRL", "ALT", "SHIFT"};

int main (int argc, char ** argv) {
	int arg=1;
	int mod=0;
	int printonly=0;
	int delay=0;
	KeySym symbol;

	if (argc < 2) usage(1,argv[0]); 
	if (!strcmp("-h",argv[arg])) usage(0,argv[0]);
	while (arg < argc-1) {
		if (!strcmp("-p",argv[arg])) { arg++; printonly=1; }
		else if (!strcmp("-d",argv[arg]) && arg < (argc-2) ) {arg++; delay=atoi(argv[arg++]); }
		else if (!strcmp("-h",argv[arg])) usage(0,argv[0]);
		else usage(1, argv[0]); 
	}
	// TODO: COMMAND LINE ARGS
	//  -p -print key-hexcode instead of sending it.
	//  -d <int> --delay <int> 
	// sleep(pre-delay);

	//
	if (!strncasecmp(argv[arg],"ALT-",4)) { mod=3; memmove(argv[arg],argv[arg]+4,strlen(argv[arg])-3); }
	if (!strncasecmp(argv[arg],"META-",4)) { mod=1; memmove(argv[arg],argv[arg]+5,strlen(argv[arg])-4); }
	if (!strncasecmp(argv[arg],"CTRL-",5)) { mod=2; memmove(argv[arg],argv[arg]+5,strlen(argv[arg])-4); }
	if (!strncasecmp(argv[arg],"SHIFT-",6)) { mod=4; memmove(argv[arg],argv[arg]+6,strlen(argv[arg])-5); }

	if (xkey_init(NULL)) { 
		fprintf(stderr,"Error: failed to connect to X server.\n");
		exit (1);
	}
	symbol=XStringToKeysym(argv[arg]);
	if (symbol==NoSymbol) {
		fprintf(stderr," Error: invalid Key-token. -- use '-h' for help.\n");
		exit(1);
	}
	if (delay) sleep(delay);
	if (printonly) printf(" key '%s' -> %lx [mod:%s]\n",argv[arg],symbol,modifier[mod]);
	else xkey_send(symbol,mod);

	return(0);
}

