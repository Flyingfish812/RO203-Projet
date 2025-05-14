import json

def sort_inner_dicts(data):
    if isinstance(data, dict):
        # 只对最内层 dict 的值是 dict 的键做排序
        return {
            k: dict(sorted(v.items())) if isinstance(v, dict) else v
            for k, v in data.items()
        }
    return data

def format_and_sort_json(input_path, output_path, indent=4):
    with open(input_path, 'r', encoding='utf-8') as infile:
        data = json.load(infile)

    # 对最内层字典排序
    sorted_data = sort_inner_dicts(data)

    with open(output_path, 'w', encoding='utf-8') as outfile:
        json.dump(sorted_data, outfile, indent=indent, ensure_ascii=False)

format_and_sort_json('solution.json', 'solution_check.json')
