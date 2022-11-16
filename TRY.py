with open('/Users/zivseker/Desktop/Projects/bio-project/xml-files/060783.xml', 'r') as result_handle:
    my_xml = result_handle.read()
    create_view_index = my_xml.find("CREATE_VIEW")
    if create_view_index != -1:
        my_xml = my_xml[:create_view_index] + my_xml[create_view_index + 11:]


with open('/Users/zivseker/Desktop/Projects/bio-project/xml-files/060783.xml', 'w') as f:
    f.write(my_xml)

