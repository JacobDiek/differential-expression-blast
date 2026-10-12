import os
import xml.etree.ElementTree as ET


def split_blast_xml(input_dir, output_dir, split_tag="Iteration", skip_already_split=False):
    os.makedirs(output_dir, exist_ok=True)
    input_xml_dir = os.listdir(input_dir)

    # Names of XML files that already have split output, so reruns only split new downloads
    already_split = set()
    if skip_already_split:
        already_split = {name.split(f"_{split_tag}_")[0] for name in os.listdir(output_dir)}

    count = 0
    failed_files = []

    for file in input_xml_dir:
        if not file.endswith(".xml"):
            continue

        file_path = os.path.join(input_dir, file)
        base_name = os.path.splitext(file)[0]
        if base_name in already_split:
            continue

        try:
            context = ET.iterparse(file_path, events=("start", "end"))
            _, root = next(context)  # get root element (e.g. <BlastOutput>)

            for event, elem in context:
                if event == "end" and elem.tag == split_tag:
                    count += 1
                    # Wrap the split element in a minimal valid XML doc
                    new_root = ET.Element(root.tag)
                    new_root.append(elem)
                    new_tree = ET.ElementTree(new_root)

                    out_path = os.path.join(output_dir, f"{base_name}_{split_tag}_{count}.xml")
                    new_tree.write(out_path, encoding="utf-8", xml_declaration=True)

                    # free memory — clear the element from the original tree
                    elem.clear()

        except ET.ParseError as e:
            print(f"Skipping {file} — parse error: {e}")
            failed_files.append(file)
            continue

    print(f"Split into {count} files in {output_dir}")
    if failed_files:
        print(f"Failed files ({len(failed_files)}): {failed_files}")
