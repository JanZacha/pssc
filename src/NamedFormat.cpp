#include "NamedFormat.h"

size_t ReplaceStringInPlace(string& subject, const string& search, const string& replace, char escape_char)
{
    /**
     * Replace all occurences of search in subject with replace except when they are preceeded by escape_char.
     * Return number of occurences.
     * Credit to: http://stackoverflow.com/questions/5343190/how-do-i-replace-all-instances-of-of-a-string-with-another-string
    */
    size_t pos = 0;
    size_t num = 0;
    while ((pos = subject.find(search, pos)) != string::npos)
    {
        if (!pos || subject[pos-1] != '%') {
            subject.replace(pos, search.length(), replace);
        }
        pos += replace.length();
        ++num;
    }
    return num;
}