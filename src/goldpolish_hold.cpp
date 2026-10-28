#include <fstream>
#include <iostream>

int
main(int argc, char** argv)
{
  if (argc < 2) {
    std::cerr << "Missing pipe path\n";
    return -1;
  }
  std::ofstream of(argv[1]);
  of << 1 << std::endl; // newline ends the message (see read_pipe_lines)
  return 0;
}
