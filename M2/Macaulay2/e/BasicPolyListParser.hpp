/**
* This class implements parsing of polynomials from a string or file
* as well as Msolve format, and writing polynomials in Msolve's input format.
*/
#ifndef M2_BASICPOLYLISTPARSER_HPP
#define M2_BASICPOLYLISTPARSER_HPP

#include <iosfwd>
#include <string>
#include <vector>

#include "BasicPolyList.hpp"

/**
 * Parses polynomials from string in the Msolve format. Msolve's format
 * includes headers for variables, a characteristic, and a list of polynomials.
 * See the 
 * <a href="https://msolve.lip6.fr/downloads/msolve-tutorial.pdf">msolve docs</a>
 * for more information information.
 *
 * \throws parsing_error
 */
BasicPolyList parseMsolveFromString(std::string contents); 
/**
* Reads the contents of the file at `filename` to a string and then calls 
* parseMsolveFromString.
 * \throws parsing_error
*/
BasicPolyList parseMsolveFile(std::string filename);

/**
 * \throws parsing_error
*/
BasicPolyList parseBasicPolyListFromString(std::string contents, std::vector<std::string> varnames);

/**
 * Generic variable names `x0, x1, ..., x(nvars-1)`, for use when writing msolve
 * input.  msolve only cares that the names are valid identifiers and that their
 * order matches the order of the variables, so there is no reason to reproduce
 * the names of the variables of the ring the polynomials came from -- in
 * general those are not valid msolve identifiers anyway, e.g. `x_(0,0)`.
 */
std::vector<std::string> msolveVarNames(int nvars);

/**
 * Writes `Fs` in msolve's *input* format, the inverse of parseMsolveFromString:
 * a line of comma separated variable names, a line with the characteristic,
 * then the polynomials, one per line, separated by commas.
 *
 * Coefficients are written as stored, so the caller is responsible for reducing
 * them into the range [0, characteristic).
 */
void writeMsolveFormat(std::ostream& o,
                       const BasicPolyList& Fs,
                       const std::vector<std::string>& varnames,
                       long characteristic);

/**
 * As writeMsolveFormat, but to the file named `filename`.
 * \throws exc::engine_error if the file cannot be opened or written
 */
void writeMsolveFile(const std::string& filename,
                     const BasicPolyList& Fs,
                     const std::vector<std::string>& varnames,
                     long characteristic);

#endif
// Local Variables:
// indent-tabs-mode: nil
// End:
