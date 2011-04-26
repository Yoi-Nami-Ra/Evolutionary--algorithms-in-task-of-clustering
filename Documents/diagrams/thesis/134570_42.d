format 75

activitynodecanvas 128042 activitynode_ref 128042 // initial_node
  xyz 27 34 2000
end
activityactioncanvas 128170 activityaction_ref 128042 // activity action Konwersja
  show_infonote default drawing_language default show_stereotype_properties default
  show_opaque_action_definition default
  xyzwh 120 13 2000 100 60
end
activityactioncanvas 128298 activityaction_ref 128170 // activity action Obliczanie odleglosci
  show_infonote default drawing_language default show_stereotype_properties default
  show_opaque_action_definition default
  xyzwh 303 13 2000 100 60
end
activityactioncanvas 128426 activityaction_ref 128298 // activity action Algorytm
  show_infonote default drawing_language default show_stereotype_properties default
  show_opaque_action_definition default
  xyzwh 478 13 2000 100 60
end
activitynodecanvas 128554 activitynode_ref 128170 // flow_final
  xyz 641 34 2000
end
flowcanvas 128682 flow_ref 128042 // <flow>
  
  from ref 128042 z 2001 to ref 128170
  show_infonote default drawing_language default show_stereotype_properties default write_horizontally default
end
flowcanvas 128810 flow_ref 128170 // <flow>
  
  from ref 128170 z 2001 to ref 128298
  show_infonote default drawing_language default show_stereotype_properties default write_horizontally default
end
flowcanvas 128938 flow_ref 128298 // <flow>
  
  from ref 128298 z 2001 to ref 128426
  show_infonote default drawing_language default show_stereotype_properties default write_horizontally default
end
flowcanvas 129066 flow_ref 128426 // <flow>
  
  from ref 128426 z 2001 to ref 128554
  show_infonote default drawing_language default show_stereotype_properties default write_horizontally default
end
end
