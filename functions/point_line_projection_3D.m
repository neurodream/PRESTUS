function [projected_point] = point_line_projection_3D(point, line_start, line_end)

line = line_end - line_start;

line_start_to_point = point - line_start;
t = dot(line_start_to_point, line) / dot(line, line);
projected_point = line_start + t * line;

end