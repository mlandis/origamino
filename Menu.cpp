/*
 * Menu.cpp
 *
 *  Created on: Jan 15, 2010
 *      Author: mikee
 */

#include "Menu.h"

Menu::Menu(Settings* sp) {
	// import settings
	settingsPtr = sp;

	// display menu
	display(sp);

}

Menu::~Menu() {
}

void Menu::display(Settings* sp) {
	// display menu
	std::cout << "Origami v0.0.1" << std::endl;
}
