from automatic_entry_blast import automatic_entry_blast
from manual_entry_blast import manual_entry_blast


def main():
    selected_tool = int(
        input(
            "What tool would you like to use?\n[1] Manual Online BLAST\n[2] Automatic Online BLAST\n[3] Automatic Manual Blast\n"
        )
    )

    if selected_tool == 1:
        manual_entry_blast()

    if selected_tool == 2:
        automatic_entry_blast()


if __name__ == "__main__":
    main()
