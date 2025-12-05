/*

Copyright 1996-2006 Roeland Merks

This file is part of Tissue Simulation Toolkit.

Tissue Simulation Toolkit is free software; you can redistribute
it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

Tissue Simulation Toolkit is distributed in the hope that it will
be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Tissue Simulation Toolkit; if not, write to the Free
Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
02110-1301 USA

*/
#include <cstdio>
#include <clocale>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <sstream>
#include <vector>
#include "sticky.h"
#include "misc.h"

/** PRIVATE **/

static int ScanForNumbers(char *string) {
    int i = 0;
    int number, numpresent = FALSE;
    char tc;
    while ((tc = string[i++]) != '\0') {
        /* tc  = char to be tested, check for end-of-string */
        number = (tc >= '0' && tc <= '9'); /* number==T if number */
        if (!(number || tc == ' ' || tc == '.')) return FALSE;
        /* return FALSE if alphanumerical char or no numbers encountered */
        numpresent |= number;
    }
    if (numpresent)
        return TRUE;
    else
        return FALSE;
}


static int ReadCleanLine(FILE *file, char *string) {
    /* Reads a line, omit the remarks */
    int c = 0;
    int result;
    char dummy;
    while (((result = fscanf(file, "%c", &string[c])) != EOF)
           && (string[c] != '\n' && string[c] != '#'))
        c++;
    /* Fast Forward until end of line */
    if (string[c] == '#')
        while (fscanf(file, "%c", &dummy) != EOF)
            if (dummy == '\n') break;

    string[c] = '\0';

    return result;
}


/** PUBLIC **/

int ReadNumber(FILE *file, int *number) {
    int result = REMARK;
    char string[255];
    setlocale(LC_NUMERIC, "C");
    while (result == REMARK) {
        if ((result = ReadCleanLine(file, string)) != EOF) {
            if (ScanForNumbers(string)) {
                sscanf(string, "%d", number);
                return OK;
            } else
                result = REMARK;
        }
    }
    return EOF;
}

int ReadDouble(FILE *file, double *number) {

    int result = REMARK;
    char string[255];
    setlocale(LC_NUMERIC, "C");

    while (result == REMARK) {
        if ((result = ReadCleanLine(file, string)) != EOF) {
            if (ScanForNumbers(string)) {
                sscanf(string, "%lf", number);
                return OK;
            } else {
                result = REMARK;
            }
        }
    }
    return EOF;
}

int FileExists(FILE **fp, const char *fname, const char *ftype) {

    *fp = fopen(fname, ftype);
    if (*fp == NULL)
        return FALSE;

    if (!strncmp(ftype, "a", 1)) {
        if (ftell(*fp) > 0L) return TRUE;
        else return FALSE;
    } else return TRUE;

}

int YesNoP(const char *message) {

    char answer[100];

    fprintf(stderr, "%s (y/n) ", message);
    fflush(stderr);

    scanf("%s", answer);
    while (strcmp(answer, "y") && strcmp(answer, "n")) {
        fprintf(stderr, "\n\bPlease answer 'y' or 'n'. ");
        fflush(stderr);
        scanf("%s", answer);
    }

    if (!strcmp(answer, "y")) return TRUE;

    return FALSE;

}

char *GetFileName(const char *message, const char *ftype) {

    static char *fname = NULL;
    FILE *tmp;
    int tmpbool;
    int must_exist;
    char tempstr[100];

    if (!strncmp(ftype, "w", 1) || !strncmp(ftype, "a", 1)) {
        must_exist = FALSE;
    } else {
        if (!strncmp(ftype, "r", 1))
            must_exist = TRUE;
        else {
            fprintf(stderr, "Error in function GetFileName: flag '%s' unknown.\n", ftype);
            exit(0);
        }
    }

    if (fname == NULL)
        fname = (char *) malloc(100 * sizeof(char));

    fprintf(stderr, "%s", message);
    fflush(stderr);

    do {
        scanf("%s", fname);
    } while ((!(tmpbool = FileExists(&tmp, fname, "r")) && must_exist
              /* File to read does not exist */
              && fprintf(stderr, "File %s not found, try again\n", fname)) ||
             /* File to write does exist */
             (tmpbool && !must_exist && sprintf(tempstr, "File %s exists, overwrite? ", fname)
              && !YesNoP(tempstr)));

    fclose(tmp);
    return fname;

}

double solveQuadradic(double a, double b, double c) {
    double root = sqrt(pow(b, 2) - 4 * a * c);
    return (-b + root) / (2 * a);
}


double dist(double x1, double y1, double x2, double y2) {
    double cat1 = (x2 - x1) * (x2 - x1);
    double cat2 = (y2 - y1) * (y2 - y1);
    return sqrt(cat1 + cat2);
}

bool pointInPolygon(std::pair<double, double> const &point, std::vector<std::pair<double, double>> const &polygon) {
    int i, j, c = 0;
    auto nvert = (int) polygon.size();
    for (i = 0, j = nvert-1; i < nvert; j = i++) {
        if ( ((polygon[i].second>point.second) != (polygon[j].second>point.second)) &&
             (point.first < (polygon[j].first-polygon[i].first) * (point.second-polygon[i].second) / (polygon[j].second-polygon[i].second) + polygon[i].first) )
            c = !c;
    }
    return c;
}

template<> std::vector<std::string> stringToVector(const std::string &values, char delim) {
    std::vector<std::string> value_vec;
    std::stringstream values_ss(values);
    std::string value_str;
    while (getline(values_ss, value_str, delim)) {
        value_vec.push_back(value_str);
    }
    return value_vec;
}