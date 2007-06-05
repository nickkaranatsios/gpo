/*
 * Copyright 2007 Google Inc.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * Author: Russell L. Smith
 */

/*
 * System configuration discovery program. The usage is:
 *
 *   configurator <config.h-file-to-generate> <compiler-command-line> <THIS_DIR-variable>
 *
 * This program operates by generating a number of test source files
 * and attempting to compile them.
 */

/***************************************************************************/
/* Standard headers. Hopefully everyone has these! */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/***************************************************************************/
/* Standard system header files */

#define NUM_HEADERS 16

char *header_file[NUM_HEADERS] = {
  "stdio.h",
  "stdlib.h",
  "string.h",
  "math.h",
  "stdarg.h",
  "assert.h",
  "errno.h",
  "fcntl.h",
  "fenv.h",
  "setjmp.h",
  "signal.h",
  "time.h",
  "unistd.h",
  "sys/stat.h",
  "sys/time.h",
  "sys/types.h"};

int has_header[NUM_HEADERS];

/*
 * Stuff that goes before all #includes
 */
char *pre_include =
  "#ifndef _LARGEFILE64_SOURCE\n"
  "#define _LARGEFILE64_SOURCE\t\t// To access lseek64() on some systems\n"
  "#endif\n";

/***************************************************************************/
/* Utility functions */

/*
 * Print an error message and exit.
 */
void fatal_error(char *message) {
  printf("\n*** Configurator failed: %s.\n\n"
	  "Please fix your configuration and try again.\n\n", message);
  exit(1);
}


/*
 * Open a file, generate an error if it can't be done.
 */
FILE * xfopen(char *filename, char *mode) {
  FILE *f;
  f = fopen(filename, mode);
  if (!f) fatal_error("Can not open a file");
  return f;
}


/*
 * Return 1 if the file exists or 0 if not.
 */
int file_exists(char *filename) {
  FILE *f;
  f = fopen(filename, "rb");
  if (f) fclose(f);
  return(f != 0);
}


/*
 * Delete a file.
 */
void delete_file(char *filename) {
  unlink(filename);
  if (file_exists(filename)) {
    fatal_error("The delete_file() function does not work");
  }
}


/*
 * Run a compile command.
 */
char *compile_cmd_line = 0;
void compile(char *output, char *input) {
  char cmd[1000];
  strcpy(cmd, compile_cmd_line);
  strcat(cmd, output);
  strcat(cmd, " ");
  strcat(cmd, input);
  printf("Running: %s\n", cmd);
  if (system(cmd) == -1) fatal_error("Could not execute command");
}


/*
 * Run a program we've just compiled, return its exit code.
 */
char *run_prefix = "";
int run(char *filename) {
  int ret;
  char cmd[1000];
  strcpy(cmd, run_prefix);
  strcat(cmd, filename);
  printf("Running: %s\n", cmd);
  ret = system(cmd);
  if (ret == -1) fatal_error("Could not execute command");
  return ret;
}


/*
 * Return 1 if a source string compiles and (if 'run_it' is nonzero) runs,
 * otherwise return 0.
 */
int compiles_and_runs(char *source, int run_it) {
  int i, ret = 0;
  FILE *f = xfopen("ctest.cpp", "wt");
  fprintf(f, "%s", pre_include);
  for (i = 0; i < NUM_HEADERS; i++) {
    if (has_header[i]) fprintf(f, "#include <%s>\n", header_file[i]);
  }
  fprintf(f, "%s", source);
  fclose(f);
  delete_file("ctest.exe");
  compile("ctest.exe", "ctest.cpp");
  if (file_exists("ctest.exe")) {
    if (run_it) {
      ret = run("ctest.exe") == 0;
    } else {
      ret = 1;
    }
  }
  delete_file("ctest.cpp");
  delete_file("ctest.exe");
  return ret;
}

/***************************************************************************/
/* Tests */

void get_standard_headers(FILE *file) {
  int i;
  fprintf(file, "\n// Standard system headers\n%s", pre_include);
  for (i = 0; i < NUM_HEADERS; i++) has_header[i] = 0;
  for (i = 0; i < NUM_HEADERS; i++) {
    has_header[i] = 1;
    if (!compiles_and_runs("int main() { return 0; }\n", 0)) has_header[i] = 0;
  }
  for (i = 0; i < NUM_HEADERS; i++) {
    if (has_header[i]) fprintf(file, "#include <%s>\n", header_file[i]);
  }
}

void get_integer_typedefs(FILE *file) {
  fprintf(file, "\n// Integer types (we assume int >= 32 bits)\n");
  if (sizeof(char) != 1) fatal_error("Expecting sizeof(char) == 1");
  if (sizeof(int) < 4) fatal_error("Expecting sizeof(int) >= 4");
  fprintf(file, "typedef char int8;\ntypedef unsigned char uint8;\n");

  if (sizeof(short) != 2) fatal_error("Can not find 2 byte integer type");
  fprintf(file, "typedef short int16;\ntypedef unsigned short uint16;\n");

  if (sizeof(short) == 4) {
    fprintf(file, "typedef short int32;\ntypedef unsigned short uint32;\n");
  } else if (sizeof(int) == 4) {
    fprintf(file, "typedef int int32;\ntypedef unsigned int uint32;\n");
  } else {
    fatal_error("Can not find 4 byte integer type");
  }

  if (sizeof(long int) == 8) {
    fprintf(file, "typedef long int64;\ntypedef unsigned long uint64;\n");
  } else {
    /*
     * 64 bit types are probably accessed via language extensions, such as
     * the gcc 'long long' type. We can't use those type directly in this
     * code because we would risk compilation errors, so we'll try compiling
     * an external program.
     */
    if (compiles_and_runs("int main() { return sizeof(long long int) != 8; }\n", 1)) {
      fprintf(file, "typedef long long int int64;\ntypedef unsigned long long int uint64;\n");
    } else {
      fatal_error("Can not find 8 byte integer type");
    }
  }

  /*
   * Get off64_t
   */
  if (!compiles_and_runs("int main() { off64_t a = 1; return 0; }\n", 0)) {
    fprintf(file, "typedef int64 off64_t;\n");
  }
}

/***************************************************************************/
/* Main */

int main(int argc, char **argv) {
  FILE *file;

  if (argc < 3 || argc > 4)
    fatal_error("Configurator expects 2 or 3 arguments");
  compile_cmd_line = argv[2];
  if (argc >= 4) run_prefix = argv[3];

  printf("\n*** Configurator starting, you may see some harmless error messages ***\n\n");

  file = xfopen(argv[1], "wt");
  fprintf(file,
    "// Per-machine configuration. This file is automatically generated.\n\n"
    "#ifndef _GPO_CONFIG_H_\n"
    "#define _GPO_CONFIG_H_\n");
  get_standard_headers(file);
  get_integer_typedefs(file);

  /*
   * See if a Pentium on gcc is available.
   */
  fprintf(file, "\n// Is this a pentium on a gcc-based platform?\n");
  fprintf(file, "#define GPO_PENTIUM_AVAILABLE %d\n", compiles_and_runs(
    "int main() {\n"
    "  asm (\"mov $0,%%eax\\n mov %%eax,(%%esi)\\n cpuid\\n\" : : : \"%eax\");\n"
    "  return 0;\n"
    "}\n", 0));

  /*
   * See if lseek64() is available.
   */
  fprintf(file, "\n// Is lseek64() available?\n");
  fprintf(file, "#define GPO_LSEEK64_AVAILABLE %d\n", compiles_and_runs(
    "int main() { lseek64(0,0,0); return 0; }\n", 0));

  /*
   * See if the O_LARGEFILE flag to open() is available.
   */
  if (!compiles_and_runs("int main() { open(\"\",O_LARGEFILE); return 0; }\n", 0)) {
    fprintf(file, "\n// This flag to open() unavailable\n#define O_LARGEFILE 0\n");
  }

  /*
   * Find or create the offsetof() macro.
   */
  fprintf(file, "\n// The offsetof() macro\n#define GPO_OFFSETOF(type, member) ");
  if (compiles_and_runs("struct A { int a,b; }; int main() { return offsetof(A,b); }\n", 0)) {
    fprintf(file, "offsetof(type,member)\n");
  } else {
    if (compiles_and_runs("struct A { int a,b; }; int main() { return __builtin_offsetof(A,b); }\n", 0)) {
      fprintf(file, "__builtin_offsetof(type, member)\n");
    } else {
      fprintf(file, "((int)(&(((type*)0)->member)))\n");
    }
  }

  fprintf(file, "\n#endif\n");
  fclose(file);

  printf("\n*** Configurator succeeded ***\n\n");
  return 0;
}
