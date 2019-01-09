/* keymap.c: watches the X keymaps and does clever things with keys */

/* Copyright (C) 1999 by the Massachusetts Institute of Technology.
 *
 * Permission to use, copy, modify, and distribute this
 * software and its documentation for any purpose and without
 * fee is hereby granted, provided that the above copyright
 * notice appear in all copies and that both that copyright
 * notice and this permission notice appear in supporting
 * documentation, and that the name of M.I.T. not be used in
 * advertising or publicity pertaining to distribution of the
 * software without specific, written prior permission.
 * M.I.T. makes no representations about the suitability of
 * this software for any purpose.  It is provided "as is"
 * without express or implied warranty.
 */

#include "nawm.h"
#include <X11/extensions/XTest.h>

extern Display *dpy;

static KeySym *syms;
static XModifierKeymap *mods;
static int keysyms_per_keycode, min_keycode, max_keycode;

static char *modname[] = { "Shift", "Lock", "Control", "Meta",
			   "Mod2", "Mod3", "Mod4", "Mod5" };

void initkeymap(void)
{
  XDisplayKeycodes(dpy, &min_keycode, &max_keycode);
  syms = XGetKeyboardMapping(dpy, min_keycode, max_keycode - min_keycode + 1,
			     &keysyms_per_keycode);
  mods = XGetModifierMapping(dpy);
}

void update_keymap(void)
{
  XFree(syms);
  XFreeModifiermap(mods);
  initkeymap();
}

struct modstate {
  KeyCode kc;
  int pressed;
};

void *set_modifier_state(int setmods)
{
  int nmods, modsize, mod, key, state, wantstate;
  struct modstate *old;
  char keys[32];
  KeyCode kc;

  nmods = 0;
  modsize = 10;
  old = xmalloc(modsize * sizeof(struct modstate));
  memset(old, 0, modsize * sizeof(struct modstate));

  XQueryKeymap(dpy, keys);

  /* For each modifier... */
  for (mod = 0; mod < 8; mod++)
    {
      wantstate = setmods & (1 << mod);

      /* For each key that generates that modifier... */
      for (key = 0; key < mods->max_keypermod; key++)
	{
	  kc = mods->modifiermap[mod * mods->max_keypermod + key];
	  if (!kc)
	    break;
	  state = KEYSTATE(keys, kc);

	  /* If the key isn't in the state we want, fix it and record
	   * that we changed it.
	   */
	  if (state != wantstate)
	    {
	      if (nmods == modsize - 1)
		{
		  old = xrealloc(old, modsize * 2 * sizeof(struct modstate));
		  memset(old + modsize * sizeof(struct modstate), 0,
			 modsize * sizeof(struct modstate));
		  modsize *= 2;
		}
	      old[nmods].kc = kc;
	      old[nmods].pressed = state;
	      nmods++;

	      if (wantstate)
		XTestFakeKeyEvent(dpy, kc, True, CurrentTime);
	      else
		XTestFakeKeyEvent(dpy, kc, False, CurrentTime);
	    }

	  /* If we want the modifier to be set, we only need to look
	   * at the first key, since we'll have set it if it wasn't
	   * already set. (If we want it unset, we need to check each
	   * key and unset any that are set.)
	   */
	  if (wantstate)
	    break;
	}
    }
  return (void *)old;
}

void reset_modifier_state(void *oldstate)
{
  struct modstate *old = (struct modstate *)oldstate;
  int i;

  for (i = 0; old[i].kc; i++)
    {
      if (!old[i].pressed)
	XTestFakeKeyEvent(dpy, old[i].kc, False, CurrentTime);
      else
	XTestFakeKeyEvent(dpy, old[i].kc, True, CurrentTime);
    }
  free(old);
}

void typechar(KeySym ks)
{
  KeyCode kc;
  int mod, xmod = -1, control = 0, start;
  void *old;

  if (ks < 32)
    {
      control = ControlMask;
      ks += 96;
    }
  kc = XKeysymToKeycode(dpy, ks);
  start = (kc - min_keycode) * keysyms_per_keycode;

  for (mod = 0; mod < 2; mod++)
    {
      /* This sucks, but OpenWindows does this wrong. (It lists
       * the uppercase letter in the 0 position and doesn't list
       * the lowercase letter at all.)
       */
      if (isalpha(ks))
	{
	  if (mod == 0 && islower(ks) && tolower(syms[start]) == ks)
	    {
	      /* Looking at unshifted position, want a lowercase
	       * letter, and this is the letter we're looking for.
	       */
	      xmod = 0;
	    }
	  else if (mod == 1 && isupper(ks) &&
		   (syms[start + mod] == ks || syms[start] == ks))
	    {
	      /* Looking at shifted position, want an uppercase
	       * letter, and either this is the letter we're looking
	       * for or the unshifted one was.
	       */
	      xmod = 1;
	    }
	}
      else if (syms[start + mod] == ks)
	xmod = mod;

      if (xmod != -1)
	{
	  old = set_modifier_state(xmod + control);
	  XTestFakeKeyEvent(dpy, kc, True, CurrentTime);
	  XTestFakeKeyEvent(dpy, kc, False, CurrentTime);
	  reset_modifier_state(old);
	  break;
	}
    }
}
