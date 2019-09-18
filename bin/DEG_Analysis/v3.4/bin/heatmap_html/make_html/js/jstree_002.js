/**
 * Created by liuzuming on 15/12/10.
 *
 * 该方法用于给用户显示文件树结构,依赖文件为:
 *  <link rel="stylesheet" type="text/css" href="<bio:contextpath/>/static/css/font-awesome/css/font-awesome.css">
 *  <link rel="stylesheet" type="text/css" href="<bio:contextpath/>/static/css/css/plugins/jsTree/style.min.css">
 *  <link href="<bio:contextpath/>/static/css/css/bootstrap.min.css" rel="stylesheet">
 *      <script src="js/jquery-2.1.1.js"></script>
 *      <script src="js/bootstrap.min.js"></script>
 *      <script src="js/plugins/jsTree/jstree.min.js"></script>
 *
 *  参数获取:
 *      tree_node_id : 生成文件树的div的id
 *      file_or_folder_type: 可选参数,暂时不使用
 *      button_click_node_id: 点击事情,点击按钮后获取该节点的原始数据
 *      http_url: 文件访问链接
 *  重写方法:
 *      jstree_global._cb_original_data(original_data)
 *      参数表示的为传入的原始对象
 *  调用方法:
 *      jstree_global.init({
            tree_node_id : "tree",
            button_click_node_id:"button_ok",
            http_url:projectName+"/util/tree/getFileTree"
        });
 *
 * 项目API网站:
 *  https://www.jstree.com/api/#/?f=get_string%28key%29
 */

;
(function () {

    var jstree_global = window.jstree_global = window.jstree_global || {};

    var _tree_node_id, _file_or_folder_type, _button_click_node_id, _http_url;

    jstree_global.init = function (core) {
        _tree_node_id = core.tree_node_id;
        if (!_tree_node_id) {
            throw "生成tree的div没有初始化";
            return;
        }
        //_file_or_folder_type = core.file_or_folder_type;
        _button_click_node_id = core.button_click_node_id;
        if (!_button_click_node_id) {
            throw "点击按钮id没有初始化";
            return;
        }
        _http_url = core.http_url;
        if (!_http_url) {
            throw "请求链接没有初始化";
            return;
        }
        _start_up_jstree();
        _add_click_event();
        return this;
    };

    /**
     * 起动jstree树形结构
     * @private
     */
    function _start_up_jstree() {
        $('#' + _tree_node_id).jstree({
            'core': {
                'data': {
                    'url': _http_url,
                    "dataType": "json",
                    "data": function (node) {
                        if (!node.original) {
                            return {"path": "#"}
                        }
                        return {"path": node.original.path}
                    }
                }
            },
            "types": {
                "default": {
                    "icon": "fa fa-folder jstree-themeicon-custom"
                }
            },
            "plugins": ["types"]
        });
    };

    /**
     * 添加按钮点击事件
     * @private
     */
    function _add_click_event() {
        $('#' + _button_click_node_id).on('click', function () {
            var tree = $('#' + _tree_node_id).jstree();
            var a = tree.get_container();
            var b = a.find("li[aria-selected=true]");
            var select_node = tree.get_node(b.attr("id"));
            var original_data = select_node.original;
            jstree_global._cb_original_data(original_data);
        });
    };
})();

/**
 * 该方法为回调函数，是调用者来自定义实现的方法
 *
 * @param original_data 传入的原始数据
 * @private
 */
jstree_global._cb_original_data = function (original_data) {
    console.log(original_data)
}
