
#ifndef FILE_H
#define FILE_H
#include <string>
#include <fstream>
#include <stdio.h>
bool file_exists (const std::string& name) {
    FILE *file = fopen(name.c_str(), "r");
    if (file!=NULL) {
        fclose(file);
        return true;
    } else {
        return false;
    }
};
#endif
