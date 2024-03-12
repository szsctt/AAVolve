import sys

def main():

    current, prev, count = "", "", 0
    for line in sys.stdin:
        line = line.strip()

        # skip empty lines
        if line == "":
            continue

        # get fields
        line_count, seq = line.split('\t')

        # sum and print result
        if seq == current:
            count += int(line_count)
        else:
            if current:
                print(f"{count}\t{current}")
            prev = current
            current = seq
            count = int(line_count)

        # check sorting
        if prev and prev > seq:
            raise ValueError(f"Input not sorted: {prev} > {seq}")

    # print last result
    if current:
        print(f"{count}\t{current}")

if __name__ == "__main__":
    main() 