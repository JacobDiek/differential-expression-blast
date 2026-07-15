import os
import xml.etree.ElementTree as ET


def split_blast_xml(input_dir, output_dir, split_tag="Iteration"):
    os.makedirs(output_dir, exist_ok=True)
    input_xml_dir = os.listdir(input_dir)    

    count = 0

    for file in input_xml_dir:
        # Grab the root's tag/attributes so each split file has a valid wrapper
        context = ET.iterparse((os.path.join(input_dir, file)), events=("start", "end"))
        _, root = next(context)  # get root element (e.g. <BlastOutput>)

        for event, elem in context:
            if event == "end" and elem.tag == split_tag:
                count += 1
                # Wrap the split element in a minimal valid XML doc
                new_root = ET.Element(root.tag)
                new_root.append(elem)
                new_tree = ET.ElementTree(new_root)

                out_path = os.path.join(output_dir, f"{split_tag}_{count}.xml")
                new_tree.write(out_path, encoding="utf-8", xml_declaration=True)

                # free memory — clear the element from the original tree
                elem.clear()

        print(f"Split into {count} files in {output_dir}")

split_blast_xml(input_dir='/Users/jacobdiek/Coding/projects/differential-expression-blast/ncbi_files_results/ncbi_files_large_alignment', 
                output_dir='/Users/jacobdiek/Coding/projects/differential-expression-blast/ncbi_files_results/ncbi_files_raw_output')