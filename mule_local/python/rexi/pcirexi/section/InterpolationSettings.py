class InterpolationSettings:
    section_number: int = 1
    interp_per_section: int = 2
    interp_method: str = 'equidistant'
    same_arc_length_for_all_sections: bool = False
    sample_nodes_based_on_arc_length: bool = True

    def __init__(self, section_number, interp_per_section, interp_method, same_arc_length_for_all_sections,
                 sample_nodes_based_on_arc_length):
        self.section_number = section_number
        self.interp_method = interp_method
        self.interp_per_section = interp_per_section
        self.same_arc_length_for_all_sections = same_arc_length_for_all_sections
        self.sample_nodes_based_on_arc_length = sample_nodes_based_on_arc_length
