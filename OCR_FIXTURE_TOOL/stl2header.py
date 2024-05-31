import os

def sanitize_name(name):
    """Sanitizes the file name to be a valid C++ variable name by replacing invalid characters."""
    return name.replace('-', '_').replace(' ', '_').replace('.', '_')

def file_to_cpp_array(input_filename, array_name):
    """Reads a binary file and converts its content to a C++ unsigned char array."""
    try:
        with open(input_filename, 'rb') as file:
            content = file.read()
        # Format each byte as a hex string
        array_content = ', '.join(f'0x{byte:02x}' for byte in content)
        return f"""const unsigned char {array_name}[] = {{
    {array_content}
}};
const size_t {array_name}_Size = sizeof({array_name});
"""
    except IOError as e:
        print(f"Error processing file {input_filename}: {e}")
        return ""
    
def files_to_cpp_header(directory, output_filename):
    """Generates a C++ header file that embeds all .stl files in a directory as arrays."""
    header_content = "#ifndef OCR_FONT_STL_H\n#define OCR_FONT_STL_H\n#pragma once\n#include <cstddef>\n#include <vector>\n\nstruct STLData {\n    std::string key;\n    const unsigned char* data;\n    size_t size;\n};\n\n"
    vector_entries = []
    for filename in sorted(os.listdir(directory)):
        if filename.lower().endswith(".stl"):
            input_filename = os.path.join(directory, filename)
            variable_base = sanitize_name(os.path.splitext(filename)[0])
            array_name = f"_{variable_base}_Data"
            header_content += file_to_cpp_array(input_filename, array_name)
            # Use the sanitized base filename as the key
            vector_entries.append(f'{{"{variable_base}", {array_name}, {array_name}_Size}}')
    # Append the vector of STLData definitions to the header content
    header_content += "\nstatic const std::vector<STLData> FONT_STL = {\n    " + ",\n    ".join(vector_entries) + "\n};\n\n#endif // OCR_FONT_STL_H\n"
    # Write the combined content to the output header file
    try:
        with open(output_filename, 'w') as file:
            file.write(header_content)
        print(f"Successfully created {output_filename}")
    except IOError as e:
        print(f"Error writing to output file: {e}")
        
if __name__ == "__main__":
    # Get the directory of the current script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    directory_path = os.path.join(script_dir, 'models')
    output_header = os.path.join(script_dir, 'OCR_font_STL.h')
    files_to_cpp_header(directory_path, output_header)
