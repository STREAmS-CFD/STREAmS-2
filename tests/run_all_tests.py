import os
import sys
import configparser
import subprocess

def update_ini_file(filename, changes):
    backup = filename + ".bak"
    os.rename(filename, backup)

    config = configparser.ConfigParser()
    config.optionxform = str
    config.read(backup)

    for section, key, value in changes:
        if not config.has_section(section):
            config.add_section(section)
        config.set(section, key, value)

    with open(filename, "w") as f:
        config.write(f)

def read_last_row(file_path):
    with open(file_path, "r") as f:
        lines = [line.strip() for line in f if line.strip()]
    return [float(x) for x in lines[-1].split()]

def diff_columns(expected, actual, tol=1e-12):
    diffs = [abs(a - b) for a, b in zip(expected, actual)]
    max_diff = max(diffs)
    if all(d == 0 for d in diffs):
        status = "OK"
    elif max_diff < tol:
        status = "Warning"
    else:
        status = "Error"
    return diffs, status

def parse_config_block(block):
    items = []
    for pair in block.split(";"):
        section, rest = pair.split(",", 1)
        key, value = rest.split("=", 1)
        items.append((section.strip(), key.strip(), value.strip()))
    return items

def read_mpi_np(ini_file):
    config = configparser.ConfigParser()
    config.optionxform = str
    config.read(ini_file)

    def clean_int(val):
        return int(val.split(";")[0].strip())

    try:
        x = clean_int(config.get("mpi", "x_split"))
        y = clean_int(config.get("mpi", "y_split"))
        z = clean_int(config.get("mpi", "z_split"))
        return x * y * z
    except Exception as e:
        print(f"Errore nella lettura dei parametri mpi: {e}")
        return 1

def main():
    TESTS_FILE = "tests.ini"
    OUT_FILE = "tests.out"
    EXECUTABLE = "../../code/streams_2.exe"
    CONFIG_NAME = "singleideal.ini"
    PROGRESS_FILE = "progress.out"

    with open(TESTS_FILE, "r") as f:
        lines = [line.strip() for line in f if line.strip()]

    results = []
    i = 0
    statuses = []
    while i < len(lines):
        if lines[i].startswith("[") and "]" in lines[i]:
            test_name = lines[i][1:-1].strip()
            param_line = lines[i + 1]
            expected_line = lines[i + 2]

            print(f"Esecuzione test: {test_name}")

            os.chdir(test_name)

            changes = parse_config_block(param_line)
            update_ini_file(CONFIG_NAME, changes)

            np = read_mpi_np(CONFIG_NAME)
            print(f"Numero di processi MPI: {np}")

            subprocess.run([
                "mpirun", "--allow-run-as-root", "--oversubscribe",
                "-np", str(np), EXECUTABLE
            ], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

            actual = read_last_row(PROGRESS_FILE)
            expected = [float(x) for x in expected_line.split()]
            diffs, status = diff_columns(expected, actual)
            statuses.append(status)

            results.append(f"[{test_name}] {status}")
            for idx, d in enumerate(diffs):
                results.append(f"  col {idx:2d}: diff = {d:.3e}")

            os.chdir("..")
            i += 3
        else:
            i += 1

    with open(OUT_FILE, "w") as f:
        f.write("\n".join(results))

    print(f"Test completati. Risultati scritti in {OUT_FILE}")
    if all(s == "OK" for s in statuses):
        sys.exit(0)
    else:
        sys.exit(1)

if __name__ == "__main__":
    main()

