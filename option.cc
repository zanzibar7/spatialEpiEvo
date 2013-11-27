#include "option.h"
#include <stdio.h>
#include <string.h>
#include <iostream>
//#include <vector>
#include <list>
#include <stdlib.h>

Option::Option(void) {}
Option::Option(char x, void* parm, const char* d, const char* e, const char* p) {
		symbol = x;
		parameter = parm;
		description = d;
		error = e;
		print_form = p;
}
Option::Option(char x, void* parm, const char* d, const char* p) {
		symbol = x;
		parameter = parm;
		description = d;
		error = d;
		print_form = p;
}
char Option::Symbol(void) { return symbol; }
const char* Option::Error(void) { return error; }
const char* Option::Description(void) { return description; }

// This is a gludy way to get different types of output
std::ostream& operator<<( std::ostream& out, const Option& a ) {
	char* v = new char[20];
	switch ( a.print_form[0] ) {
		case 'd': // double
			sprintf(v,"%g",*((double*) a.parameter));
			break;
		case 'i': // int
			sprintf(v,"%d",*((int*) a.parameter));
			break;
		case 'b': // bool
			if ( *((bool *) a.parameter ) == true )
			{ sprintf(v,"true"); }
			else
			{ sprintf(v,"false"); }
			break;
		case 'c': // char
			sprintf(v,"%c",*((char*) a.parameter));
			break;
		case 's': // string
			sprintf(v,"%s",*((char**) a.parameter));
			break;
		default:
			std::cerr << a.print_form << " not printable " << std::endl;
			exit(1); // What should this be?
	}
	out << v ;
	delete v;
	return out;
}

OptionBool::OptionBool(char x, void* parm, const char* d, const char* e):
		Option(x, parm, d, e, "b") {};
OptionBool::OptionBool(char x, void* parm, const char* d):
		Option(x, parm, d, "b") {};
bool OptionBool::read(const char arg[]) {
	//std::cerr << arg << '\t' << parameter << std::endl;
	const char * s = &(arg[2]);
	if ( s[0] == '0' )
	{ *((bool *)parameter) = false; }
	else if ( s[0] == '1' )
	{ *((bool *)parameter) = true; }
	else
	{ *((bool *)parameter) = ! *((bool *)parameter); }
	return true;
}
std::ostream& OptionBool::operator<<( std::ostream& out ) {
	if ( *((bool *)parameter) == true ) 
	{ out << "true"; }
	else
	{ out << "false"; }
	return out;
}

OptionDouble::OptionDouble(char x, void* parm, const char* d, const char* e):
	Option(x, parm, d, e, "d") {};
OptionDouble::OptionDouble(char x, void* parm, const char* d):
	Option(x, parm, d, "d") {};
bool OptionDouble::read(const char* arg) {
	const char * s = &(arg[2]);
	int i = sscanf(s,"%lf",parameter);
	if ( 0 == i || EOF == i ) return false;
	else return true;
}

OptionInt::OptionInt(char x, void* parm, const char* d, const char* e):
		Option(x, parm, d, e, "i") {};
OptionInt::OptionInt(char x, void* parm, const char* d):
		Option(x, parm, d, "i") {};
bool OptionInt::read(const char* arg) {
	const char * s = &(arg[2]);
	int i = sscanf(s,"%d",parameter);
	// in the upgrade to gcc 4.3, the flag from reading an int
	// seems to have changed from %ld to %d?
	if ( 0 == i || EOF == i ) return false;
	else return true;
}
//std::ostream& OptionInt::operator<<( std::ostream& out ) {
//	out << *((int *) parameter);
//	return out;
//}

OptionChar::OptionChar(char x, void* parm, const char* d, const char* e):
	Option(x, parm, d, e, "c") {};
OptionChar::OptionChar(char x, void* parm, const char* d):
	Option(x, parm, d, "c") {};
bool OptionChar::read(const char* arg) {
	const char * s = &(arg[2]);
	int i = sscanf(s,"%c",parameter);
	if ( 0 == i || EOF == i ) return false;
	else return true;
}


OptionString::OptionString(char x, void* parm,
		const char* d, const char* e):
		Option(x, parm, d, e, "s") {};
OptionString::OptionString(char x, void* parm,
		const char* d):
		Option(x, parm, d, "s") {};
bool OptionString::read(const char* arg) {
	const char * s = &(arg[2]);
	int i = sscanf(s,"%s",parameter);
	if ( 0 == i || EOF == i ) return false;
	else return true;
}
//std::ostream& OptionString::operator<<( std::ostream& out ) {
//	out << (char *) parameter;
//	return out;
//}

void Options::Assert(bool t, const char* s) {
	if ( false == t ) {
		std::cerr << "Command line error: " << s << std::endl;
		exit(1); // What should this be?
	}
}
Options::Options(void) {
	message = (const char**)NULL;
	exit_on_verbose = true;
};

Options::Options(const char* m[]) {
	message = m;
	exit_on_verbose = true;
};

void Options::AddOption(Option* x) { list_of_options.push_back(x); }

void Options::AddIntOption(char x, void* parm, const char* d, const char* e) {
	OptionInt* o = new OptionInt(x, parm, d, e);
	list_of_options.push_back(o);
}
void Options::AddDoubleOption(char x, void* parm, const char* d, const char* e) {
	OptionDouble* o = new OptionDouble(x, parm, d, e);
	list_of_options.push_back(o);
}
void Options::AddBoolOption(char x, void* parm, const char* d, const char* e) {
	OptionBool* o = new OptionBool(x, parm, d, e);
	list_of_options.push_back(o);
}

void Options::CheckOptions(const int argc, const char** argv) {
	int mark;
	bool found;
	for (mark=1; argv[mark] && argv[mark][0]=='-';mark++) {
		found = false;
		for (std::list<Option*>::iterator j=list_of_options.begin(); j!=list_of_options.end(); ++j){
			if ( argv[mark][1] == (*j)->Symbol() ) {
				Assert((*j)->read(argv[mark]), (*j)->Error());
				found = true;
				break;
			}
		} if ( found ) continue;
		if (argv[mark][1] == 'v' || argv[mark][1] == '?') {
			for (int i=0; message && message[i]; ++i) {
				std::cerr << message[i] << std::endl;
			}
			for (std::list<Option*>::iterator j=list_of_options.begin(); j!=list_of_options.end(); ++j) {
				//std::cerr << *gji << std::endl;
				//std::cerr << (*j)->Description() << std::endl;
				//std::cerr << (*j)->Description() << *gji << std::endl;
				std::cerr << (*j)->Description() << " [" << **j << "]" << std::endl;
			}
			if (exit_on_verbose) exit(0);
			else continue;
		}
		std::cerr << "Ignored argument: " << argv[mark] << std::endl;
	}
}
void Options::Print(void) {
	for (std::list<Option*>::iterator j=list_of_options.begin(); j!=list_of_options.end(); ++j){
		std::cout
		<< "# " << **j  << (*j)->Description() << std::endl
		;
	}
}
