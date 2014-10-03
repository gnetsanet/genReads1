#ifndef STRREF_H
#define STRREF_H

class StringRef
{
private:
	char const*     begin_;
	int             size_;

public:
	int size() const { return size_; }
	char const* begin() const { return begin_; }
	char const* end() const { return begin_ + size_; }

	StringRef( char const* const begin, int const size );
};

#endif