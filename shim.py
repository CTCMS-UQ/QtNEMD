#!/usr/bin/python3
import lammps
import InputManager

def main():
    lmp = lammps.lammps()
    params = InputManager.InputManager()
    params.reset_and_update(lmp)
    print(lmp.numpy.extract_atom("x"))
    print(lmp.numpy.extract_atom("y"))
    #lmp.command("run 10 pre yes post no")
    #lmp.command("run 1000 pre no post yes")

if __name__ == "__main__":
    main()
