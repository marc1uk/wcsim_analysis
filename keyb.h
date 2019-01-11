/*
kbhit() and getch() for Linux/UNIX
Chris Giese <geezer@execpc.com>	http://my.execpc.com/~geezer
*/

#define LINUX  // or pass -DLINUX to g++
#ifdef LINUX

/*****************************************************************************/
/*  GETCH  */
/*****************************************************************************/
int getch(void);


/*****************************************************************************/
/*  KBHIT  */
/*****************************************************************************/
int kbhit();

#endif
