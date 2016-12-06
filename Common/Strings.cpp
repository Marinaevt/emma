#include "StdAfx.h"
#include "strings.h"

//---------------------------------------------------------------------------

StringVector split(const std::string& s, size_t max, const std::string& delims)
{
        StringVector splits;
        size_t num = 0,
                 start = 0,
                 pos;

        do
        {
                pos = s.find_first_of(delims, start);

                if (pos == start)
                {
                        start = pos + 1;
                }
                else if (pos == s.npos || (max && num == max))
                {
                        // copy the rest of the string
                        splits.push_back(s.substr(start));
                        break;
                }
                else
                {
                        // copy up to delimiter
                        splits.push_back(s.substr(start, pos - start));
                        start = pos + 1;
                }

                // parse up to next real data
                start = s.find_first_not_of(delims, start);
                ++num;
        }
        while (pos != s.npos);

        return splits;
}
//---------------------------------------------------------------------------

std::string trim(const std::string& s, bool left, bool right)
{
        size_t i, j;

        i = (left ? s.find_first_not_of(" \t\n") : 0);
        j = (right ? s.find_last_not_of(" \t\n") : s.length());

        return s.substr(i, j - i + 1);
}
//---------------------------------------------------------------------------

double toDouble(const std::string& s)
{
        return atof(s.c_str());
}
//---------------------------------------------------------------------------

int toInt(const std::string& s)
{
        return atoi(s.c_str());
}
//---------------------------------------------------------------------------

bool toBool(const std::string& s)
{
		if (s == "false" || s == "0" || s.empty()) {
			return false;
		}
		return true;
}
//---------------------------------------------------------------------------