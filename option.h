#ifndef _option_cc
#define _option_cc
// A collection of classes for easy command-line
// parsing.  Currently can handle boolean, int,
// and double types. And maybe string?
/*
 * To use it, first create a variable opt of class options.
 * Then, repeatedly call AddOption with appropriate arguements
 */

//#include <vector>
#include <list>
#include <iostream>

class Option {
	private:
	protected:
		char *name;
		const char *description;
		const char *error;
		const char *print_form;
		char symbol;
		void *parameter;
	public:
		Option();
		Option(char x, void* parm, const char* d, const char* p);
		Option(char x, void* parm, const char* d, const char* e, const char* p);
		virtual bool read(const char* arg)=0;
		//virtual bool read(const char* arg)=false;
		
		//virtual std::ostream& operator<<  ( std::ostream& );
		//virtual ostream& operator<< (ostream&, const Option& );
		char Symbol(void);
		const char* Error(void);
		const char* Description(void);
		friend std::ostream& operator<< (std::ostream&, const Option&);
};

class Options {
	private:
		//std::vector<Option*> v; // deleted 2005-11-27
		std::list<Option*> list_of_options;
		void Assert(bool t, const char* s);
		const char** message;
	public:
		bool exit_on_verbose;
		Options(void);
		Options(const char**);
		void AddOption(Option* x);

		void AddIntOption(char x, void* parm, const char* d, const char* e);
		void AddDoubleOption(char x, void* parm, const char* d, const char* e);
		void AddBoolOption(char x, void* parm, const char* d, const char* e);

		void CheckOptions(const int argc, const char** argv);
		void Print(void);
};

class OptionBool : public Option {
	public:
		OptionBool(char x, void* parm, const char* d, const char* e);
		OptionBool(char x, void* parm, const char* d);
		bool read(const char* arg);
		std::ostream& operator<<( std::ostream& out );
};

class OptionDouble : public Option {
	public:
		OptionDouble(char x, void* parm, const char* d, const char* e);
		OptionDouble(char x, void* parm, const char* d);
		bool read(const char* arg);
		std::ostream& operator<<( std::ostream& out );
};

class OptionInt : public Option {
	public:
		OptionInt(char x, void* parm, const char* d, const char* e);
		OptionInt(char x, void* parm, const char* d);
		bool read(const char* arg);
		std::ostream& operator<<( std::ostream& out );
};

class OptionChar : public Option {
	public:
		OptionChar(char x, void* parm, const char* d, const char* e);
		OptionChar(char x, void* parm, const char* d);
		bool read(const char* arg);
		std::ostream& operator<<( std::ostream& out );
};

class OptionString : public Option {
	public:
		OptionString(char x, void* parm, const char* d, const char* e);
		OptionString(char x, void* parm, const char* d);
		bool read(const char* arg);
		std::ostream& operator<<( std::ostream& out );
};

// ostream& operator<< (ostream& __o, const Option& __r)

#endif
