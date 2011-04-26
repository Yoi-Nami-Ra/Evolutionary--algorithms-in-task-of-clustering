format 75

activitynodecanvas 128024 activitynode_ref 128024 // initial_node
  xyz 13 12 2000
end
activitynodecanvas 128042 activitynode_ref 141610 // decision
  xyz 10 75 2000
end
activitynodecanvas 128408 activitynode_ref 128280 // flow_final
  xyz 216 83 2000
end
note 128426 "Czy blok zewnetrzny lub watek zewn trzny ?"
  xyzwh 43 10 2006 141 53
activitynodecanvas 128810 activitynode_ref 141738 // decision
  xyz 11 241 2000
end
activityactioncanvas 129066 activityaction_ref 141610 // activity action zaladuj rekord wiersza
  show_infonote default drawing_language default show_stereotype_properties default
  show_opaque_action_definition default
  xyzwh 158 190 2000 100 60
end
note 129322 "Czy rekord pierwszej kolumny ?"
  xyzwh 71 107 2006 121 49
activitynodecanvas 129578 activitynode_ref 141866 // decision
  xyz 195 283 2000
end
activitynodecanvas 129962 activitynode_ref 141994 // decision
  xyz 195 382 2000
end
note 130218 "Czy blok skrajny ?"
  xyzwh 16 302 2000 111 35
note 130474 "Czy blok pierwszego wiersza ?"
  xyzwh 19 361 2000 79 61
activityactioncanvas 130730 activityaction_ref 141738 // activity action zaladuj rekord kolumny
  show_infonote default drawing_language default show_stereotype_properties default
  show_opaque_action_definition default
  xyzwh 265 372 2000 100 60
end
activitynodecanvas 131114 activitynode_ref 142250 // join
  horizontal xyz 303 255 2000
end
activityactioncanvas 132266 activityaction_ref 141866 // activity action Oblicz odleglosc
  show_infonote default drawing_language default show_stereotype_properties default
  show_opaque_action_definition default
  xyzwh 265 62 2000 100 60
end
note 132906 "Synchronizuj watki"
  xyzwh 204 128 2000 89 53
flowcanvas 128170 flow_ref 142506 // <flow>
  
  from ref 128024 z 2001 to ref 128042
  show_infonote default drawing_language default show_stereotype_properties default write_horizontally default
end
flowcanvas 128298 flow_ref 142634 // tak
  
  from ref 128042 z 2001 label "tak" xyz 117 77 3000 to ref 128408
  show_infonote default drawing_language default show_stereotype_properties default write_horizontally default
end
flowcanvas 128938 flow_ref 142762 // <flow>
  
  from ref 128042 z 2001 to ref 128810
  show_infonote default drawing_language default show_stereotype_properties default write_horizontally default
end
flowcanvas 129194 flow_ref 142890 // tak
  
  from ref 128810 z 2001 label "tak" xyz 89 230 3000 to ref 129066
  show_infonote default drawing_language default show_stereotype_properties default write_horizontally default
end
flowcanvas 129706 flow_ref 143018 // <flow>
  
  from ref 128810 z 2001 to ref 129578
  show_infonote default drawing_language default show_stereotype_properties default write_horizontally default
end
flowcanvas 129834 flow_ref 143146 // <flow>
  
  from ref 129066 z 2001 to ref 129578
  show_infonote default drawing_language default show_stereotype_properties default write_horizontally default
end
flowcanvas 130090 flow_ref 143274 // nie
  
  from ref 129578 z 2001 label "nie" xyz 200 342.5 3000 to ref 129962
  show_infonote default drawing_language default show_stereotype_properties default write_horizontally default
end
flowcanvas 130858 flow_ref 143402 // tak
  
  from ref 129962 z 2001 label "tak" xyz 225.5 383 3000 to ref 130730
  show_infonote default drawing_language default show_stereotype_properties default write_horizontally default
end
flowcanvas 131242 flow_ref 143530 // <flow>
  decenter_begin 243
  
  from ref 129578 z 2001 to ref 131114
  show_infonote default drawing_language default show_stereotype_properties default write_horizontally default
end
flowcanvas 131754 flow_ref 143786 // <flow>
  
  from ref 130730 z 2001 to ref 131114
  show_infonote default drawing_language default show_stereotype_properties default write_horizontally default
end
flowcanvas 131882 flow_ref 143914 // <flow>
  
  from ref 129962 z 2001 to ref 131114
  show_infonote default drawing_language default show_stereotype_properties default write_horizontally default
end
flowcanvas 132394 flow_ref 144042 // <flow>
  
  from ref 131114 z 2001 to ref 132266
  show_infonote default drawing_language default show_stereotype_properties default write_horizontally default
end
flowcanvas 132522 flow_ref 144170 // <flow>
  
  from ref 132266 z 2001 to ref 128408
  show_infonote default drawing_language default show_stereotype_properties default write_horizontally default
end
line 128682 -_-_
  from ref 128426 z 2007 to ref 128042
line 129450 -_-_
  from ref 128810 z 2007 to ref 129322
line 130346 -_-_ decenter_begin 513
  from ref 129578 z 2001 to ref 130218
line 130602 -_-_ decenter_begin 654
  from ref 130474 z 2001 to ref 129962
line 133034 -_-_
  from ref 132906 z 2001 to ref 131114
end
