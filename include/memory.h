#ifndef MEMORY_H_
#define MEMORY_H_

#include <algorithm>
#include <cstdio>
#include <cstring>

using namespace std;

inline int parseLine(char* line) {
  int i = strlen(line);
  while (*line < '0' || *line > '9') line++;
  line[i - 3] = '\0';
  i = atoi(line);
  return i;
}

inline int getValueOfPhysicalMemoryUsage(){
  FILE* file = fopen("/proc/self/status", "r");
  int result = -1;
  char line[128];

  while (fgets(line, 128, file) != NULL){
    if (strncmp(line, "VmRSS:", 6) == 0){
      result = parseLine(line);
      break;
    }
  }
  fclose(file);
  return result;
}

#endif
