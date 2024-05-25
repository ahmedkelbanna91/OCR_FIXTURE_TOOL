import os

def sanitize_name(name):
    """Sanitizes the file name to be a valid C++ variable name."""
    return name.replace('-', '_').replace(' ', '_')

def files_to_cpp_header(directory, output_filename):
    header_content = "#pragma once\n\n"
    
    # Loop through all .stl files in the directory
    for filename in sorted(os.listdir(directory)):
        if filename.lower().endswith(".stl"):
            input_filename = os.path.join(directory, filename)
            # Create a valid C++ variable name for the array
            array_name_base = sanitize_name(os.path.splitext(filename)[0])
            array_name = f"_{array_name_base}"
            array_size_name = f"{array_name}_Size"
            
            try:
                # Read the binary content of the input file
                with open(input_filename, 'rb') as file:
                    content = file.read()

                # Convert the binary content to a C++ unsigned char array
                array_content = ', '.join(f'0x{byte:02x}' for byte in content)

                # Append the header content for the current file
                header_content += f"""const unsigned char {array_name}[] = {{
    {array_content}
}};
const size_t {array_size_name} = sizeof({array_name});

"""
            except IOError as e:
                print(f"Error processing file {filename}: {e}")

    # Write the combined header file
    try:
        with open(output_filename, 'w') as file:
            file.write(header_content)
        print(f"Successfully created {output_filename}")
    except IOError as e:
        print(f"Error writing to output file: {e}")

# Example usage - Replace 'your_directory_path' with the path to your STL files
directory_path = 'C:/Users/kline/Desktop/DEV/OCR_STL/x64/Release/new-ocr-a'
output_header = 'C:/Users/kline/Desktop/DEV/OCR_STL/OCR_STL/new-ocr-a.h'
files_to_cpp_header(directory_path, output_header)
