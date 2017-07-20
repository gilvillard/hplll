
// ************************************
// From http://www.latticechallenge.org
// ************************************

#ifndef TOOLS_H
#define TOOLS_H

template<typename T>
std::istream& operator>>(const std::string& in, T& v)
{
    static std::istringstream* iss = 0;
    if (iss) delete iss;
    iss = new std::istringstream(in);
    *iss>>v;
    return *iss;
}

#define PARSE_MAIN_ARGS \
    std::ostringstream message; \
    message << argv[0];  \
    for (int i=0; i<argc; ++i)

#define MATCH_MAIN_ARGID(arg,dest) \
    if (i==0)  {\
	message << " [" << arg << " " << dest << "]"; \
    } else { \
	if (std::string(argv[i]) == std::string(arg))  {\
	    std::string(argv[i+1]) >> dest; \
	    i++; continue; \
	} \
    }

#define DETECT_MAIN_ARGID(arg,dest,value) \
    if (i==0)  {\
	message << "[" << arg << "]"; \
    } else { \
	if (std::string(argv[i]) == std::string(arg))  {\
	    dest = value; \
	    continue; \
	} \
    }

#define SYNTAX() \
    if (i==0) continue; \
    std::cerr << "Syntax: " << message.str() << std::endl; \
    exit(1);


#endif
