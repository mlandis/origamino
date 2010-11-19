/*
 * Menu.h
 *
 *  Created on: Jan 15, 2010
 *      Author: mikee
 */

#ifndef MENU_H_
#define MENU_H_

#include <iostream>

#include "Settings.h"

class Settings;

class Menu {

public:
	Menu(Settings*);
	virtual ~Menu();

private:
	void display(Settings*);
	Settings* settingsPtr;

};

#endif /* MENU_H_ */
