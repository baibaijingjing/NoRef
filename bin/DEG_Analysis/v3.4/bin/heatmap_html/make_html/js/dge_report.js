/*window.DgeReport.zhuliucheng.asdfs*/
;
(function(window,$,appUtil, undefined){

    var PageStyle,   //页面样式问题
        _isBlank,   //是否为空结题报告
        _rootPath=$("#rootPath").val(),//结题报告根目录
        _projectPath,//项目相对路劲
        _projectName, //项目名
        _drawMapPath="/share/bioCloud/cloud/temp/",//画图的临时目录
        _projectId, //结题报告编号
        _typeproject,//项目版本
        projectCode, //项目编号
        dataContent = {},//相对路径配置文件
        Util, //工具模块
        _maxTabCount = 6, //数据挖掘页面最多显示的标签个数
        geneName = "",//表达基因选择的对比的样品
        isCytoscape = false,//是否制作了蛋白图
        DataMining,//数据挖掘模块
        ProteinIntMapping,//蛋白互作图

        folder="",//无参功能注释数据导入，选择的文件夹
        gene_file="",//无参Unigene数据导入，选择的文件
        dgeselect="",//获得差异基因挖掘的差异分组

        GFF,
        Genome_unigene,
        Preparation,

        oneOrTwo="multiple",//判断是单端数据还是双端数据
        typeselect,
        parameter="you",

        error=false,

        folderSin="",fq1Sin="",folderMul="",fq1Mul="",fq2Mul="",

        Process,//主流程模块
        _hasProp = Object.prototype.hasOwnProperty,
        _extends = function (child, parent) {
        for (var key in parent) {
            if (_hasProp.call(parent, key)) {
                child[key] = parent[key];
            }
        }
        return child;
        },
        _bind = function (fn, me) {
            return (function () {
                fn.apply(me, arguments);
            });
        },
        Dge=window.DgeReport=window.DgeReport||{};    //对象挂载在windows下
        var jstree_global = window.jstree_global = window.jstree_global || {};

    //数据初始化
    Dge.init = function (initData) {
        if (!(_isBlank = initData.isBlank)) { // 默认为该页面没有结题报告
            _isBlank = false;
        }
        if (!(_rootPath = initData.rootPath) && !_isBlank) {
            throw "结题报告更目录没有初始化，请重试！";
            return;
        }
        if (!(_projectName = initData.projectName)) {
            if (_projectName == "") {
            } else {
                throw "项目名没有初始化，请重试！project name no init !";
                return;
            }
        }
        if (!_isBlank && !(_drawMapPath = initData.drawMapPath)) {
            throw "绘图临时目录没有初始化，请重试！drawMapPath no init ";
            return;
        }
        if (initData.projectId) {
            _projectId = initData.projectId;
        }
        if (initData.maxTabCount) {
            _maxTabCount = initData.maxTabCount;
        }
        return this;
    };
    //开始加载页面数据和样式
    Dge.run = function () {
        for (var o in PageStyle.initialCode) {
            var value = PageStyle.initialCode[o];
            if (typeof value === "function") {
                value();
            }
        }
        Process.initReferenceSpecies();
        Process.nextSteps();
        Dge.FDR_onclick();
        Dge.FC_onclick();
        $(".clearfix a[href='#finish']").hide();
        $(".clearfix a[href='#cancel']").hide();
        if (!_isBlank) {
            $.ajax({type : "POST",
                async:false,
                cache:false,
                url :"/report/dge/reportCommon",
                dataType:"html",
                data:{"projectId":_projectId},
                success : function (data){
                    if(data.message =="success"){
                        $("#685782A07D2D425DB87BEEADD9352700").append(data);
                    }
                },
                error:function(err){
                    console.log(err);
                }
            });
            try {
                Process.init();
            } catch (e) {
                console.log(e);
            }
            DataMining.init();
            DataMining.DegGroupSelect.openPop();
            DataMining.DegGroupSelect.submit();
            Process.DegGroupingSelect.pop();
            Process.initReferenceSpeciesFromDetail();
        }
        return this;
    };
    //FDR FC页面初始化触发事件
    Dge.FDR_onclick =function(){
        $("input[name='FDR']").click(function(){
            if($("input[name='FDR']:checked").val()=="Customer"){
                $("input[name='FDR'][type='text']").css("display","inline-block");
            }else{
                $("input[name='FDR'][type='text']").css("display","none");
            }
        });
    }
    Dge.FC_onclick =function(){
        $("input[name='FC']").click(function(){
            if($("input[name='FC']:checked").val()=="Customer"){
                $("input[name='FC'][type='text']").css("display","inline-block");
            }else{
                $("input[name='FC'][type='text']").css("display","none");
            }
        });
    }
    //公共模块
    Util = (function () {
        var Util = new Object();
        _extends(Util, appUtil);

        Util.proxy = function (func) {
            return _bind(function () {
                return func.apply(this, arguments);
            }, this);
        };
        Array.prototype.indexOf = function(elt , from)
        {
            var len = this.length >>> 0;
            var from = Number(arguments[1]) || 0;
            from = (from < 0)
                ? Math.ceil(from)
                : Math.floor(from);
            if (from < 0)
                from += len;
            for (; from < len; from++)
            {
                if (from in this &&
                    this[from] === elt)
                    return from;
            }
            return -1;
        };
        return Util;
    })();
    //主流程模块
    Process = (function () {
        var Process = new Object();
        //选择有参或者无参
        Process.show=function(_this) {
            $("button[name=selectType]").addClass("btn-outline");
            $(_this).removeClass("btn-outline");
            if($(_this).val()=="shen"){
                parameter="you";
                $("#shen").show();
                $("#no_reference").hide();
                $("#7C57F3B43FC04B6A9CE3DD91FDCDB534").show();
            }else if($(_this).val()=="no_reference"){
                parameter="wu";
                $("#shen").hide();
                $("#no_reference").show();
                $("#7C57F3B43FC04B6A9CE3DD91FDCDB534").hide();
            }
        }
        _extends(Process, Util);
        //初始化参考物种
        Process.initReferenceSpecies = function(){
            var refSp = $("#7C8CEE4D2E264AD2B3D67601B0DA1A5A"),
                data = appUtil.postJSON(_projectName+"/refReport/RNARefReport/referenceSpecies"),
                gene = $("#11523B85157341A4A31F118243F568B2"),cookie=document.cookie,
                html=""
                ,geneHtml="";
            if(data ==null){
                return;
            }
            if(data ==null){
                return;
            }
            var num = 0;
            var i=0;
            if(cookie.indexOf("locale=zh_CN")<0){//英文
                i=1;
            }else {//中文
                i=0;
            }
            if(data ==null){
                return;
            }
            if(data ==null){
                return;
            }
            var num = 0;
            for(var item in data ){
                if(item.indexOf("_1")<0){
                    var str ="";
                    for(var itemtow in data){
                        if(itemtow == item || itemtow.indexOf(item)>=0 ){
                            str += data[itemtow].split("\t")[2]+",";
                        }
                    }
                    if(num == 0 ){
                        html = "<div class='alert alert-info checked' onclick='showMessage(this)' style='background-color: #1AB394'><a href='javascript:void(0);' " +
                        "class='alert-link' style='color:#FFF'>"+data[item].split("\t")[i]+"</a>"+
                        "<input type='hidden' value='"+str.substring(0,str.length-1)+"'/>"+
                        "</div>";
                        geneHtml="<div class='alert alert-primary alert-dismissable js-1'> "+
                        "<input type='radio' checked style='float: left;margin-left: 10px;margin-top: 4px;'  value='"+data[itemtow].split("\t")[i+2]+"'/>"+
                        "<a style='margin-left: 15px' href='javascript:void(0);'  class='alert-link'>"+data[item].split("\t")[i+2]+"</a>"+
                        "</div>";
                        num++;
                    }else{
                        html += "<div class='alert alert-info checked' onclick='showMessage(this)' ><a href='javascript:void(0);' " +
                        "class='alert-link'>"+data[item].split("\t")[i]+"</a>"+
                        "<input type='hidden' value='"+str.substring(0,str.length-1)+"'/>"+
                        "</div>";
                    }
                }
            }
            gene.html(geneHtml);
            refSp.html(html);
            Process.referenceSpecies();
            Process.referenceSpecies = data;
        }
        //初始化默认值
        var detailCfg = {
            "Project_name":"default",
            "Customer_info":"default",
            "Project_id":"default",
            "Project_key":"default",
            "Customer_name": "Zhang",
            "First_time": "2015/02/13",
            "Second_time": "2015/02/13",
            "Third_time": "2015/02/13",
            "Contract_data": "8M",
            "Q30": "85%",
            "Lib_type":"fr-unstranded",
            "Mismatch":"2",
            "Insert_size":"40",
            "Memory":"15G",
            "Queue_type":"middle.q",
            "Com":"",
            "fold":"2",
            "FDR":"0.01"
        };
        var annotation = {
            "Bacteria": {
                "nr": "/share/nas2/database/ncbi/Nt_Nr_division/nr_BCT",
                "Kegg": "/share/nas2/database/kegg/kegg_divide/Bacteria.fa",
                "Cog": "/share/nas2/database/cog/myva",
                //"Kog": "/share/nas2/database/kog/kyva",
                //"Pfam": "/share/nas2/database/pfam/27.0/Pfam-A.hmm",
                "Swissprot": "/share/nas2/database/uniprot/knowledgebase/current/TrEMBL_Divide/uniprot_sprot_bacteria.fasta"
            },
            "Plants": {
                "nr": "/share/nas2/database/ncbi/Nt_Nr_division/nr_PLN",
                "Kegg": "/share/nas2/database/kegg/kegg_divide/Plants.fa",
                "Cog": "/share/nas2/database/cog/myva",
                //"Kog": "/share/nas2/database/kog/kyva",
                //"Pfam": "/share/nas2/database/pfam/27.0/Pfam-A.hmm",
                "Swissprot": "/share/nas2/database/uniprot/knowledgebase/current/TrEMBL_Divide/uniprot_sprot_plants.fasta"
            },
            "Invertebrates": {
                "nr": "/share/nas2/database/ncbi/Nt_Nr_division/nr_INV",
                "Kegg": "/share/nas2/database/kegg/kegg_divide/Animals.fa",
                "Cog": "/share/nas2/database/cog/myva",
                //"Kog": "/share/nas2/database/kog/kyva",
                //"Pfam": "/share/nas2/database/pfam/27.0/Pfam-A.hmm",
                "Swissprot": "/share/nas2/database/uniprot/knowledgebase/current/TrEMBL_Divide/uniprot_sprot_invertebrates.fasta"
            },
            "Vertebrates": {
                "nr": "/share/nas2/database/ncbi/Nt_Nr_division/nr_VRT",
                "Kegg": "/share/nas2/database/kegg/kegg_divide/Animals.fa",
                "Cog": "/share/nas2/database/cog/myva",
                //"Kog": "/share/nas2/database/kog/kyva",
                //"Pfam": "/share/nas2/database/pfam/27.0/Pfam-A.hmm",
                "Swissprot": "/share/nas2/database/uniprot/knowledgebase/current/TrEMBL_Divide/uniprot_sprot_vertebrates.fasta"
            },
            "Primates": {
                "nr": "/share/nas2/database/ncbi/Nt_Nr_division/nr_PRI",
                "Kegg": "/share/nas2/database/kegg/kegg_divide/Animals.fa",
                "Cog": "/share/nas2/database/cog/myva",
                //"Kog": "/share/nas2/database/kog/kyva",
                //"Pfam": "/share/nas2/database/pfam/27.0/Pfam-A.hmm",
                "Swissprot": "/share/nas2/database/uniprot/knowledgebase/current/TrEMBL_Divide/uniprot_sprot_mammals.fasta"
            },
            "Rodents": {
                "nr": "/share/nas2/database/ncbi/Nt_Nr_division/nr_ROD",
                "Kegg": "/share/nas2/database/kegg/kegg_divide/Animals.fa",
                "Cog": "/share/nas2/database/cog/myva",
                //"Kog": "/share/nas2/database/kog/kyva",
                //"Pfam": "/share/nas2/database/pfam/27.0/Pfam-A.hmm",
                "Swissprot": "/share/nas2/database/uniprot/knowledgebase/current/TrEMBL_Divide/uniprot_sprot_rodents.fasta"
            },
            "Mammals": {
                "nr": "/share/nas2/database/ncbi/Nt_Nr_division/nr_MAM",
                "Kegg": "/share/nas2/database/kegg/kegg_divide/Animals.fa",
                "Cog": "/share/nas2/database/cog/myva",
                //"Kog": "/share/nas2/database/kog/kyva",
                //"Pfam": "/share/nas2/database/pfam/27.0/Pfam-A.hmm",
                "Swissprot": "/share/nas2/database/uniprot/knowledgebase/current/TrEMBL_Divide/uniprot_sprot_mammals.fasta"
            },
            "Fungi": {
                "nr": "/share/nas2/database/ncbi/nr",
                "Kegg": "/share/nas2/database/kegg/kegg_divide/Fungi.fa",
                "Cog": "/share/nas2/database/cog/myva",
                //"Kog": "/share/nas2/database/kog/kyva",
                //"Pfam": "/share/nas2/database/pfam/27.0/Pfam-A.hmm",
                "Swissprot": "/share/nas2/database/uniprot/knowledgebase/current/TrEMBL_Divide/uniprot_sprot_fungi.fasta"
            },
            "Total": {
                "nr": "/share/nas2/database/ncbi/nr",
                "Kegg": "/share/nas2/database/kegg/genes.pep",
                "Cog": "/share/nas2/database/cog/myva",
                //"Kog": "/share/nas2/database/kog/kyva",
                //"Pfam": "/share/nas2/database/pfam/27.0/Pfam-A.hmm",
                "Swissprot": "/share/nas2/database/uniprot/knowledgebase/current/complete/uniprot_sprot.fasta"
            }
        };

        var dataCfg;
        var referenceSpecies,
            FCdefined, FDRdefined, SepDefined,
            ComDefined, para_K_covDefined, annoDatabseType, annoDatabseTypeNew,fileSize="";
        Process.init = function () {
            _initCfgFile.apply(this, null);
            OriginalData.init(); // 用来初始化原始数据
            if("Com" in detailCfg){
                ComDefined = _assignmentTovar(detailCfg.Com);
            }
            if("Sep" in detailCfg){
                SepDefined = _assignmentTovar(detailCfg.Sep);
            }
            FCdefined = _assignmentTovar(detailCfg.fold);
            FDRdefined = _assignmentTovar(detailCfg.FDR);
            this.proxy(_initPage());
        };


        //参考物种确定
        Process.submitReferSpecies=function(){
            if(parameter=="you"){
                /*var refSp=$("#7C8CEE4D2E264AD2B3D67601B0DA1A5A").find("div.checked").find("a").text(),json = Process.referenceSpecies,
                    gene=$("#11523B85157341A4A31F118243F568B2").find("input[type='radio']").is(":checked");

                if(json && refSp){
                    for(var item in json){
                        if(item == refSp || item == refSp+"_1"){
                            var detail = json[item].split("\t");
                            detailCfg.Genome_unigene = detail[5].split(":")[1];
                            detailCfg.Preparation = detail[4].split(":")[1];
                            detailCfg.GFF = detail[6].split(":")[1];
                        }
                    }
                }
                if(!gene){
                    appUtil.popbox ("请选择参考物种及组装版本！");
                    return false;
                }else{
                    return true;
                }*/

                var json = Process.referenceSpecies,gen=$("#11523B85157341A4A31F118243F568B2").find("input[type='radio']");
                var gene,g;
                for(var i=0;i<gen.length;i++){
                    if(($(gen[i])).is(":checked")){
                        console.info(($(gen[i])));
                        g=($(gen[i]));
                        gene=($(gen[i])).parent(".alert-primary").find("a").text();
                    }
                }
                if(json ){
                    for(var item in json){
                        console.info(json[item].indexOf(gene));
                        if(json[item].indexOf(gene)>=0  || json[item].indexOf(gene+"_1")>=0 ){
                            var detail = json[item].split("\t");
                            detailCfg.Genome_unigene = detail[5].split(":")[1];
                            detailCfg.Preparation = detail[4].split(":")[1];
                            detailCfg.GFF = detail[6].split(":")[1];
                        }
                    }
                }
                if(!g){
                    appUtil.popbox("请选择参考物种及组装版本！");
                    return false;
                }else{
                    return true;
                }
            }else if(parameter=="wu"){
                var s1=$("#button_ok1").val();//文件
                var s2=$("#userName1").val();//文件夹
                if(s1==""){
                    appUtil.popbox ("Unigene数据不能为空，请选择数据！");
                    return false;
                }else if(s1.contains(".fa")){
                    return true;
                }else {
                    appUtil.popbox ("您选择的文件格式不是.fa格式的，请重新选择！");
                    return false;
                }
                if(s2==""){
                    appUtil.popbox ("功能注释数据不能为空，请选择数据！");
                    return false;
                }
            }
        }
        //参考物种选择
        Process.referenceSpecies = function(){
            $("#7C8CEE4D2E264AD2B3D67601B0DA1A5A").change(function(){
                var gene = $("#11523B85157341A4A31F118243F568B2"),json = Process.referenceSpecies,
                    refSp = $("#7C8CEE4D2E264AD2B3D67601B0DA1A5A").val(),html="";

                if(json && refSp){
                    for(var item in json){
                        if(item == refSp || item == refSp+"_1"){
                            var str = json[item].split("\t");
                            html += "<option value="+str[2]+">"+str[2]+"</option>";
                        }
                    }
                    gene.html(html);
                }
            });
        }
        //参考物种初始化回显  返回选择,没有添加none
        Process.initReferenceSpeciesFromDetail= function(){
            var int = 0,cookie=document.cookie;
            var j=0;
            if(cookie.indexOf("locale=zh_CN")<0){//英文
                j=1;
            }else {//中文
                j=0;
            };
            if(detailCfg.Genome_unigene && detailCfg.Preparation && detailCfg.GFF  && Process.referenceSpecies){
                for(var items in Process.referenceSpecies){
                    var refSpecies = Process.referenceSpecies[items];
                    if(refSpecies.indexOf(detailCfg.Genome_unigene)>0 && refSpecies.indexOf(detailCfg.Preparation)>0
                        && refSpecies.indexOf(detailCfg.GFF)>0){
                        //Process.referenceSpecies();
                        $("#7C8CEE4D2E264AD2B3D67601B0DA1A5A").val(refSpecies.split("\t")[1]);

                        var aRef=$("#7C8CEE4D2E264AD2B3D67601B0DA1A5A").find("a");
                        var aGene=$("#11523B85157341A4A31F118243F568B2")/*.find("a")*/;
                        for(var i=0;i< aRef.length;i++){
                            var initRefSpecies= refSpecies.split("\t")[j];//中文"("+refSpecies.split("\t")[0]+")"
                            if(initRefSpecies== aRef[i].text){
                                /*选中的样式*/
                                $(aRef[i]).parent("div.alert-info").css('background-color','#1AB394');
                                $(aRef[i]).parent("div.alert-info").addClass("checked");
                                $(aRef[i]).parent("div.alert-info").find("a").css('color','#fff');
                                /*未选中的样式*/
                                $(aRef[i]).parent("div.alert-info").siblings().removeClass("checked");
                                $(aRef[i]).parent("div.alert-info").siblings().css('background-color','#fff');
                                $(aRef[i]).parent("div.alert-info").siblings("div.alert-info").find("a").css({'color':'#676A6C'});

                            }
                        }
                        $(aGene).html("<div class='alert alert-primary alert-dismissable js-1 '>"+
                        "<input type='radio' checked style='float: left;margin-left: 10px;margin-top: 4px;'  value='"+refSpecies.split("\t")[2]+"'>"+
                        "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<a href='javascript:void(0);' id='assemble1' style='color:#fff' class='alert-link'>"+refSpecies.split("\t")[2] +"</a>"+
                        "</div>");
                        Process.refSpeciesData= refSpecies;
                        int ++;
                    }
                }
            }
        }
        //用来初始化两个配置文件
        function _initCfgFile() {
            var data = this.postJSON("/report/dge/getCfg", {projectId: _projectId});
            if (data.message) {
                detailCfg = Process.detailCfgNext = data.detail;
                dataCfg = Process.dataCfgNext = data.data;
                typeselect=data.typeSelect == null ? "" : data.typeSelect;
                if(typeselect=="single"){
                    $("#single").attr("checked",true);
                    oneOrTwo="single";
                }else if(typeselect=="multiple"){
                    $("#multiple").attr("checked",true);
                    oneOrTwo="multiple";
                }
                if(data.selected=="no"){
                    parameter="wu";
                    $("#wu").removeClass("btn-outline");
                    $("#you").addClass("btn-outline");
                    $("#shen").hide();
                    $("#no_reference").show();
                    $("#7C57F3B43FC04B6A9CE3DD91FDCDB534").hide();
                }else if(data.selected=="yes"){
                    parameter="you";
                    GFF=detailCfg.GFF;
                    Genome_unigene= detailCfg.Genome_unigene;
                    Preparation=detailCfg.Preparation;

                    $("#wu").addClass("btn-outline");
                    $("#you").removeClass("btn-outline");
                    $("#shen").show();
                    $("#no_reference").hide();
                    $("#7C57F3B43FC04B6A9CE3DD91FDCDB534").show();
                }
                //打开原始数据，显示之前有项目入口的值
                var sampleList=[],item;
                for(item in data.data){
                    var obj={};
                    obj.name=item;
                    if (obj.name!="Qphred"){
                        sampleList.push(obj);
                    }
                }
                Process.samples=sampleList.sort(function (a, b) {//样品名称排序
                    var a1 = parseInt(a.name.substring(1));
                    var b1 = parseInt(b.name.substring(1));
                    return a1-b1;
                });
            }else{
                appUtil.popbox (data.message);
            }
        };
        //将变量clone，方便后来对比，参数时候改变
        function _assignmentTovar(param2) {
            if (!param2) {
                throw "配置文件错误，无法操作";
                return;
            } else {
                return param2;
            }
        };
        //初始化页面
        function _initPage() {
            pageShow.FDRSelect(detailCfg.FDR);
            pageShow.FCSelect(detailCfg.fold);
            if("Project_key" in detailCfg){
                pageShow.speciesName(detailCfg.Project_key);
            }else{
                pageShow.speciesName("Demo");
            }
            var o = pageShow.annoDatabseType();
            if (o) {
                annoDatabseType = o;
                this.annoDatabseType = o;
                pageShow.showAnnoDatabseType(o);
            } else {
                throw new Error("没有该物种的基因注释");
            }
        };
        var pageShow = (function () {
            pageShow = new Object();
            //初始化物种名称
            pageShow.speciesName = function(value){
                $("#titInput").val(value == null ? 'Demo' : value);
            }
            //FDR选项页面变化
            pageShow.FDRSelect = function (value) {
                var options = $("#FDR").find("option");
                var flag = false;
                for (var i = 0; i < options.length; i++) {
                    if ($(options[i]).val() == value) {
                        $(options[i]).attr("selected", true);
                        flag = true;
                        continue;
                    }
                }
                if (!flag) {
                    //$("#FDR").append("<option value=\"" + value + "\" selected=\"true\">" + value + "</option>");
                }
            };
            //FDR选项页面变化
            pageShow.FCSelect = function (value) {
                var asdf = $("#FC").find("option[value=" + value + "]").attr("value");
                if (!asdf) {
                    /*$("#FC").append("<option value=\"" + value + "\">" + value + "</option>");
                    $("#FC").find("option[value=" + value + "]").attr("selected", true);*/
                } else {
                    $("#FC").find("option[value=" + value + "]").attr("selected", true);
                }
            };
            //注释物种
            pageShow.annoDatabseType = function () {
                for (var o in annotation) {
                    var dg = detailCfg, an = annotation[o];

                    if (dg.nr=== an.nr &&
                        dg.Kegg === an.Kegg &&
                        dg.Swissprot === an.Swissprot
                        && dg.Cog === an.Cog ) {
                        return o;
                    }
                }
                return null;
            };
            //基因注释 页面变化
            pageShow.showAnnoDatabseType = function (value) {
                $("#annoDatabseType").find("option").removeAttr("selected");
                $("#annoDatabseType").find("option[value=" + value + "]").attr("selected", true);
                Process.annoDatabseType = value;
            };
            return pageShow;

        })();
        //FDR和差异筛选倍数阈值输入值校验
        Process.FDRfocusInput = function (_this) {
            var ref = /^0\.[0-9]{1,5}?$/;/*/^\d+(?=\.{0,1}\d+$|$)/;*/
            if (!ref.test($(_this).val())) {
                $(_this).siblings(".FDRtitInput_content").show();
            } else {
                $(_this).siblings(".FDRtitInput_content").hide();
            }
        }
        //FC输入值校验
        Process.FCfocusInput = function (_this) {
            var ref = new RegExp('^[0-9]+$', 'g');
            if (!ref.test($(_this).val())) {
                $(_this).siblings(".FCtitInput_content").show();
            } else {
                $(_this).siblings(".FCtitInput_content").hide();
            }
        }
        //原始数据选择模块代码
        var OriginalData = (function () {
            var OriginalData = OriginalData || {};
            //var OriginalData=new Object();
            _extends(OriginalData, Util);
            var els = {
                "popClick": "#originalDataImportPop", //原始数据弹出框触发按钮
                "popHtml": "#originalDataImportPopHtml",  //原始数据弹出框
                "liHtml": "#originalDataLi", //新加入的样品信息html代码
                "sampleId": "#a31c8e2c1feb4d509806c3a0bfbe2b9e", //样品信息id
                "fold": "#99f384b299944d2c894196d861b686b3", //数据文件夹
                "file1": "#5019bac9fd2a4f66b4b8496584789313", //手动输入文件1
                "file2": "#d5cc2aad50a54138abbf8084a84a1476" //手动输入文件2
            };
            //物种名称校验
            OriginalData.focusInput = function (_this) {
                if($("#titInput").val()==""){
                    appUtil.popbox ("物种名称不能为空，请重新输入!");
                    $("#titInput").css("background-color","#FBE3E4");
                    $("#titInput").css({"border":"1px solid FBC2C4"});
                    return false;
                }else{
                    var ref = new RegExp('^[A-Za-z_0-9]{1,15}$');
                    if (!ref.test($("#titInput").val())) {
                        appUtil.popbox ("物种名称输入有误，请重新输入!");
                        $("#titInput").css("background-color","#FBE3E4");
                        $("#titInput").css({"border":"1px solid FBC2C4"});
                        return false;
                    } else {
                        $("#titInput").css("background-color","#ffffff");
                        $("#titInput").css({"border":"1px solid 1AB394"});
                        return true;
                    }
                }
            }
            //项目名称校验
            OriginalData.projectInput = function () {
                //var ss;
                if($("#project_ID").val()==""){
                    appUtil.popbox ("项目名称不能为空，请重新输入!");
                    $("#project_ID").css("background-color","#FBE3E4");
                    $("#project_ID").css({"border":"1px solid FBC2C4"});
                    //ss= "no";
                    return false;
                }else{
                    var ref = new RegExp('^[\u4E00-\u9FA5A-Za-z0-9_]{1,50}$');
                    if (!ref.test($("#project_ID").val())) {
                        appUtil.popbox ("项目名称输入有误，请重新输入!");
                        $("#project_ID").css("background-color","#FBE3E4");
                        $("#project_ID").css({"border":"1px solid FBC2C4"});
                        //ss= "no";
                        return false;
                    } else {
                        $("#project_ID").css("background-color","#ffffff");
                        $("#project_ID").css({"border":"1px solid 1AB394"});
                        //ss= "yes";
                        return true;
                    }
                }
            }

            //主流程第三步next
            Process.annoSpeciesNext = function(){
                if(OriginalData.projectInput()){//校验
                    //添加差异分组
                    DegGroupingSelect.degGoupSelect();
                    return true;
                }else{
                    return false;
                }
            }
            //主流程第四步next
            Process.differentialNext = function(){
                //FDR,FC校验
                if(Process.FDRCheck()&& DegGroupingSelect.submit()){

                    return true;
                }else{
                    return false;
                }
            }

            Process.FDRCheck =function(){
                FDRdefined=$("input[name='FDR']:checked").val();
                if(FDRdefined=="Customer"){
                    var FDR=/^0\.[0-9]{1,5}?$/;/*/^\d+(?=\.{0,1}\d+$|$)/;*/
                    FDRdefined=$("input[name='FDR'][type='text']").val();
                    if( !FDR.test(FDRdefined)){
                        appUtil.popbox ("您输入的'FDR的值'不正确,请您调整后再提交任务.(注:FDR的范围为(0,0.99999])");
                        return false;
                    }
                }
                FCdefined=$("input[name='FC']:checked").val();
                if(FCdefined=="Customer"){
                    //var FC= /^\d+(?=\.{0,1}\d+$|$)/;
                    var reg = new RegExp("^[0-9]+$");
                    FCdefined=$("input[name='FC'][type='text']").val();
                    if( !reg.test(FCdefined)){
                        appUtil.popbox ("您输入的'差异筛选倍数阈值'不正确,请您调整后再提交任务。(注:差异筛选倍数阈值的范围是 (1,1000])");
                        return false;
                    }else if (FCdefined <= 1 || FCdefined > 1000) {
                        appUtil.popbox ("请输入正确的FC数值,范围是(1,1000]");
                        return "";
                    }
                }
                return true;
            }

            //判断data.cfg中的文件名字和dataContent中相同,则返回文件名字
            function _getNameFromCfg(name, no) {
                var nameCon = "";
                for (var i in dataContent) {
                    if (i == name) {
                        if (no == 2) {
                            nameCon = dataContent[i].filename2;
                        } else {
                            nameCon = dataContent[i].filename1;
                        }
                    }
                }
                return nameCon;
            }

            //原始数据初始化  名称从dataContent中取,优先从dataContent取
            OriginalData.init = function () {
                var cfg = dataCfg;
                if (!cfg) {
                    throw "配置文件出错，无法操作";
                    return;
                } else {
                    this.originalDataSin = this.originalDataSin || {items: []};//单端数据的原始数据
                    this.originalDataMul = this.originalDataMul || {items: []};//双端数据的原始数据
                    var arrSin = [];
                    var arrMul = [];
                    for (var o in cfg) {
                        if(o == "Qphred"){
                            continue;
                        }
                        var obj = cfg[o];
                        var json = {};
                        if(typeselect=="single"||oneOrTwo=="single"){
                            var path = obj["fq1"];
                            //if(_getNameFromCfg(o)=="")
                            JSON.stringify(dataContent).length == 2 ? json["fileName1"] = path.substring(path.lastIndexOf("\/") + 1) : json["fileName1"] = _getNameFromCfg(o);
                            //json["fileName1"] = path.substring(path.lastIndexOf("\/") + 1) ;
                            json["fq1"] = obj["fq1"];
                            json.name = o;
                            json.ok = true;
                            json.id = this.getUID();
                            arrSin[arrSin.length] = json;
                        }else if(typeselect=="multiple"||oneOrTwo=="multiple"){
                            var path = obj["fq1"];
                            JSON.stringify(dataContent).length == 2 ? json["fileName1"] = path.substring(path.lastIndexOf("\/") + 1) : json["fileName1"] = _getNameFromCfg(o);
                            //json["fileName1"] = path.substring(path.lastIndexOf("\/") + 1);
                            json["fq1"] = obj["fq1"];
                            var path = obj["fq2"];
                            JSON.stringify(dataContent).length == 2 ? json["fileName2"] = path.substring(path.lastIndexOf("\/") + 1) : json["fileName2"] = _getNameFromCfg(o);
                            //json["fileName2"] = path.substring(path.lastIndexOf("\/") + 1);
                            json["fq2"] = obj["fq2"];
                            json.name = o;
                            json.ok = true;
                            json.id = this.getUID();
                            arrMul[arrMul.length] = json;
                        }
                    }
                    this.originalDataSin.items.push.apply(this.originalDataSin.items, arrSin);
                    this.originalDataMul.items.push.apply(this.originalDataMul.items, arrMul);
                    this.originalDataSin.oldSize = this.originalDataSin.size = this.originalDataSin.items.length || 0;
                    this.originalDataMul.oldSize = this.originalDataMul.size = this.originalDataMul.items.length || 0;
                }
                if(typeselect=="single"){
                    //初始化原始数据导入页面样品
                    Process.originalSampleData=this.originalDataSin;
                    $("#originalDataLi1").tmpl(this.originalDataSin).appendTo(".SampleInformation");
                }else if(typeselect=="multiple"){
                    Process.originalSampleData=this.originalDataMul;
                    $("#originalDataLi").tmpl(this.originalDataMul).appendTo(".SampleInformation");
                }
                OriginalData.focusInput();
                OriginalData.projectInput();
            };
            //原始数据的数据回显
            OriginalData.pop = function () {
                    if(!_isBlank){
                        if(!this.originalDataSin||!this.originalDataMul){
                            appUtil.popbox ("没有初始化项目参数!");
                        }
                    }
                    if(typeselect=="single"){
                        $("#originalDataLi1").tmpl(this.originalDataSin).appendTo(".SampleInformation");
                        $("#single").attr("checked",true);
                        $("#data_fq2").hide();
                    }else if(typeselect=="multiple"){
                        $("#originalDataLi").tmpl(this.originalDataMul).appendTo(".SampleInformation");
                        $("#multiple").attr("checked",true);
                        $("#data_fq2").show();
                    }
                    if($("#single").is(":checked")==true){
                        oneOrTwo="single";
                    }else if($("#multiple").is(":checked")==true){
                        oneOrTwo="multiple";
                    }
            };
            //单端模式和双端模式的切换
            OriginalData.oneSelectTwo=function(ss){
                if(ss==1){
                    if($("#single").is(":checked")==true){
                        oneOrTwo="single";
                        $(".SampleInformation").html("");
                        $("#originalDataLi1").tmpl(this.originalDataSin).appendTo(".SampleInformation");
                        $("#data_fq2").hide();
                        $("#button_ok").val(folderSin);
                        $("#userName").val(fq1Sin);
                    }
                }else if(ss==2){
                    if($("#multiple").is(":checked")==true){
                        oneOrTwo="multiple";
                        $(".SampleInformation").html("");
                        $("#originalDataLi").tmpl(this.originalDataMul).appendTo(".SampleInformation");
                        $("#data_fq2").show();
                        $("#button_ok").val(folderMul);
                        $("#userName").val(fq1Mul);
                        $("#fq2Id").val(fq2Mul);
                    }
                }
            }
            //将移除的id放入该数组中   id是items中的id
            OriginalData.removeOriginalDataLi = function (id) {
                this.removeOriginalDataList = this.removeOriginalDataList || [];
                this.removeOriginalDataList.push(id);
                $("#a31c8e2c1feb4d509806c3a0bfbe2b9e").html("");
                if($("#single").is(":checked")){
                    $("#originalDataLi1").tmpl(this.originalDataSin).appendTo(".SampleInformation");
                }else if($("#multiple").is(":checked")){
                    $("#originalDataLi").tmpl(this.originalDataMul).appendTo(".SampleInformation");
                }
                OriginalData.submit();
                if(oneOrTwo=="single"){
                    var items=Process.OriginalData.originalDataSin.items;
                    for(var item in items){
                        if(id==items[item].id){
                            items.splice(item,1);
                            Process.OriginalData.originalDataSin.items.splice(item,1);
                        }
                    }
                }else if(oneOrTwo=="multiple"){
                    var items=Process.OriginalData.originalDataMul.items;
                    for(var item in items){
                        if(id==items[item].id){
                            items.splice(item,1);
                            Process.OriginalData.originalDataMul.items.splice(item,1);
                        }
                    }
                }
            };
            //判断某一li是否被删除
            OriginalData.containsId = function (id) {
                if (!this.removeOriginalDataList) {
                    return false;
                }
                var l = this.removeOriginalDataList.length;
                for (var i = 0; i < l; i++) {
                    if (id == this.removeOriginalDataList[i]) {
                        return true;
                    }
                }
                return false;
            };
            //原始数据提交按钮
            OriginalData.submit = function () {
                    var num;
                    if(oneOrTwo==""){
                        if(typeselect=="single"){
                            num=1;
                        }else if(typeselect=="multiple"){
                            num=2;
                        }
                    }else if(oneOrTwo=="single"){
                        num=1;
                    }else if(oneOrTwo=="multiple"){
                        num=2;
                    }
                    if (num==1) {
                        var subDataCfg = {};
                        var arr = this.originalDataSin.items;
                        var length = arr.length;
                        for (var i = 0; i < length; i++) {
                            if (!this.containsId(arr[i]["id"])) {
                                subDataCfg[arr[i]["name"]] = {};
                                arr[i].ok = true;
                                subDataCfg[arr[i]["name"]]["fq1"] = arr[i]["fq1"];
                            } else {
                                delete dataContent[arr[i]["name"]];
                                arr.splice(i, 1);
                                length--;
                                i--;
                            }
                        }
                        dataCfg = subDataCfg;
                        _modifyDeg();
                        this.removeOriginalDataList = [];
                        this.originalDataSin.oldSize = this.originalDataSin.size;
                        OriginalData.oriClosePop();
                    }
                    if (num==2) {
                        var subDataCfg = {};
                        var arr = this.originalDataMul.items;
                        var length = arr.length;
                        for (var i = 0; i < length; i++) {
                            if (!this.containsId(arr[i]["id"])) {
                                subDataCfg[arr[i]["name"]] = {};
                                arr[i].ok = true;
                                subDataCfg[arr[i]["name"]]["fq1"] = arr[i]["fq1"];
                                subDataCfg[arr[i]["name"]]["fq2"] = arr[i]["fq2"];
                            } else {
                                delete dataContent[arr[i]["name"]];
                                arr.splice(i, 1);
                                length--;
                                i--;
                            }
                        }
                        dataCfg = subDataCfg;
                        _modifyDeg();
                        this.removeOriginalDataList = [];
                        this.originalDataMul.oldSize = this.originalDataMul.size;
                        OriginalData.oriClosePop();
                    }

                    if (dataCfg.length == 0) {
                        throw "差异分组的配置文件没有初始化，请初始化在试？";
                        return;
                    } else {
                        var gUID = this.getUID();
                        this.gUID = gUID;  //用来存放确定按钮的id
                        var jqueryTmplData = {"items": [], "gUID": gUID, "results": []};
                        var oriData;
                        if(oneOrTwo==""){
                            if(typeselect=="single"){
                                oriData=Process.OriginalData.originalDataSin.items;
                            }else if(typeselect=="multiple"){
                                oriData=Process.OriginalData.originalDataMul.items;
                            }
                        }else if(oneOrTwo=="single"){
                            oriData=Process.OriginalData.originalDataSin.items;
                        }else if(oneOrTwo=="multiple"){
                            oriData=Process.OriginalData.originalDataMul.items;
                        }
                        var str="";
                        for (var i = 0; i < oriData.length; i++) {
                            var item = oriData[i];
                            for (var name in item) {
                                if (name == 'name') {
                                    str+=item[name]+",";
                                }
                            }
                        }
                        detailCfg.Kmean=str.substring(0,str.length-1);
                    }
            };
            //判断输入的样本名称是否相同
            function _checkInputSampleName(){
                var sampleDom = $("#a31c8e2c1feb4d509806c3a0bfbe2b9e").find("div");
                var sampleArr = "",bool =false ;
                $(sampleDom).each(function(){
                    if(sampleArr.indexOf($(this).find("input").val())>=0){
                        bool =  true;
                    }
                    sampleArr += $(this).find("input").val()+";";
                });
                return bool;
            }
            //原始数据修改后，有关差异分组进行关联修改
            function _modifyDeg() {
                var dataArr = [];
                for (var name in dataCfg) {
                    dataArr[dataArr.length] = name;
                }
                if (SepDefined) {
                    var ssss1 = SepDefined.split("!");
                    var sepStr = "";
                    for (var m = 0; m < ssss1.length; m++) {
                        if (_removeUserlessitem(dataArr, ssss1[m])) {
                            (sepStr == "") ? (sepStr += ssss1[m]) : (sepStr += "!" + ssss1[m]);
                        }
                    }
                    SepDefined = sepStr;
                }
                if (ComDefined) {
                    var ssss1 = ComDefined.split("!");
                    var comStr = "";
                    for (var m = 0; m < ssss1.length; m++) {
                        if (_removeUserlessitem(dataArr, ssss1[m])) {
                            (comStr == "") ? (comStr += ssss1[m]) : (comStr += "!" + ssss1[m]);
                        }
                    }
                    ComDefined = comStr;
                }
            };
            //old选择文件弹出框
            OriginalData.openProjectById = function () {
                var treeObj = $.fn.zTree.getZTreeObj("tree");
                var type = $("#tree").attr("name");
                if (type == "123") {
                    var fileId = treeObj.getSelectedNodes()[0].id;
                    window.location.href ="/report/report/intoByFile?fileId=" + fileId;
                }else if (type == "5") {  //无参功能注释数据导入，选择文件夹
                    _oriSelect.apply(this, [type, treeObj]);
                }else if (type == "6") {  //无参Unigene数据导入，选择文件
                    gene_file=treeObj.getSelectedNodes()[0].path;
                } else{
                    _oriSelect.apply(this, [type, treeObj]);
                }
                removeFile();
            };
            //判断添加的样品路径名称是否在样品信息中  在返回true  不再返回false
            function _isOriSelect(json) {
                var items;
                if(oneOrTwo=="single"){
                    items=Process.OriginalData.originalDataSin.items;
                }else if(oneOrTwo=="multiple"){
                    items=Process.OriginalData.originalDataMul.items;
                }

                for (var i in items) {
                    var item = items[i];
                    for (var m in item) {
                        if (m == 'fq1') {
                            if (item[m] == json) {
                                return true;
                            }
                        }
                    }
                }
            }

            //取名字的前面部分
            OriginalData._namePrefix = function (string) {
                if (string) {
                    var reg = new RegExp("[^A-Za-z0-9]", 'g'),
                        reg2 = new RegExp("[1-9*]fq", "g");
                    string = string.replace(reg, "").replace(reg2, "");
                    return string;
                }
            }
            //原始数据导入
            OriginalData.treeShow = function(id,type){
                //如果fa文件为空的时候，不能选择文件夹
                if(id=="folder_id1"){
                    var str=$("#button_ok1").val().trim();
                    if(str==""){
                        appUtil.popbox ("您导入的Unigene数据为空，请导入后选择功能注释数据！！");
                        return;
                    }
                }
                if( $("#"+id).children().length > 0){
                    if("none" == $("#"+id).css("display")){
                        $("#"+id).show();
                        $("#"+id).parent().parent().parent().siblings().find(".tree_show_id").hide();
                    }else{
                        $("#"+id).html("");
                        $("#"+id).hide();
                    }
                } else{
                    //初始化原始数据导入
                    var jqueryTmplData = {"tree_user":id+"tree_user_weiz","tree_submit":id+"tree_submit_filed"};
                    $("#tree_tmpl_html").tmpl(jqueryTmplData).prependTo($("#"+id));
                    jstree_global.init({
                        tree_node_id : id+"tree_user_weiz",
                        button_click_node_id:id+"tree_submit_filed",
                        http_url:"/report/dge/getFileTree"
                    });
                    jstree_global._cb_original_data = function (original_data) {
                        var jsonPath=appUtil.postJSON("/report/util/getUserPath");
                        //判断提交方式
                        if("folder" == type){
                            if($("#single").is(":checked")){
                                if(original_data.fileType=="folder"){
                                    var datas = appUtil.postJSON("/report/dge/getFileTree/",{path:original_data.path,filetype:"folder1"});
                                    if(JSON.stringify(datas)=="{}"){
                                        appUtil.popbox ("您导入的文件夹内，没有符合单端的数据，请重新导入！！");
                                        return;
                                    }else{
                                        folderSin=original_data.text;
                                        $("#button_ok").val(original_data.path.replace(jsonPath.userRootPath,"/根目录"));
                                        $(".SampleInformation").html("");
                                        _oriSelect(datas,"1");
                                    }
                                }else{
                                    appUtil.popbox ("您导入数据不是文件夹，请重新导入！！");
                                    return;
                                }


                            }else if($("#multiple").is(":checked")){
                                if(original_data.fileType=="folder"){
                                    var datas = appUtil.postJSON("/report/dge/getFileTree/",{path:original_data.path,filetype:"folder"});
                                    if(datas==""){
                                        appUtil.popbox ("您导入的文件夹内，没有符合双端的数据，请重新导入！！");
                                        return;
                                    }else{
                                        folderMul=original_data.text;
                                        $("#button_ok").val(original_data.path.replace(jsonPath.userRootPath,"/根目录"));
                                        $(".SampleInformation").html("");
                                        _oriSelect(datas,"1");
                                    }
                                }else{
                                    appUtil.popbox ("您导入数据不是文件夹，请重新导入！！");
                                    return;
                                }
                            }
                        }else{
                            if ("fastq1" == type) {
                                if($("#single").is(":checked")){
                                    if(original_data.fileType=="fastq"){
                                        if((original_data.text).contains("1.fq")){
                                            var  num=parseInt(1000*Math.random());
                                            var arr = [];
                                            var json = {};
                                            json.id=original_data.id;
                                            json.name="S0"+num;
                                            json["fileName1"]=original_data.text;
                                            json["fq1"] = original_data.path;
                                            json.path1=original_data.path;
                                            json.ok = true;
                                            arr[arr.length] = json;
                                            _getFile_from(arr);
                                            $("#userName").val(original_data.path.replace(jsonPath.userRootPath,"/根目录"));
                                            fq1Sin=original_data.text;
                                        }else{
                                            appUtil.popbox ("您导入的可能是双端测序数据，请将数据类型设定为双端！！");
                                            return;
                                        }
                                    }else{
                                        appUtil.popbox ("您导入数据不是fq文件，请重新导入！！");
                                        return;
                                    }

                                }else if($("#multiple").is(":checked")){
                                    if(original_data.fileType=="fastq"){
                                        var datas = appUtil.postJSON("/report/dge/getFileTree/", {path: original_data.path, filetype: "fastq"});
                                        console.log("获取单个文件" + JSON.stringify(datas));
                                        if(datas==""){
                                            appUtil.popbox ("您导入数据不符合双端数据的格式，请重新导入！！");
                                            return;
                                        }else{
                                            $("#userName").val(original_data.path.replace(jsonPath.userRootPath,"/根目录"));
                                            fq1Mul=original_data.text;
                                            _oriSelect( datas,type, original_data.path);
                                        }
                                    }else{
                                        appUtil.popbox ("您导入数据不是fq文件，请重新导入！！");
                                        return;
                                    }

                                }
                                OriginalData.cansel_tree();
                            }else if("fastq2" == type){
                                if(original_data.fileType=="fastq"){
                                    var datas = appUtil.postJSON("/report/dge/getFileTree/", {path: original_data.path, filetype: "fastq"});
                                    console.log("获取单个文件" + JSON.stringify(datas));
                                    if(datas==""){
                                        appUtil.popbox ("您导入数据不符合双端数据的格式，请重新导入！！");
                                        return;
                                    }else{
                                        $("#fq2Id").val(original_data.path.replace(jsonPath.userRootPath,"/根目录"));
                                        fq2Mul=original_data.text;
                                        _oriSelect(datas, type, original_data.path);
                                        OriginalData.cansel_tree();
                                    }
                                }else{
                                    appUtil.popbox ("您导入数据不是fq文件，请重新导入！！");
                                    return;
                                }

                            }else if("fasta" == type){
                                if(original_data.fileType=="fasta"){
                                    $("#button_ok1").val(original_data.path.replace(jsonPath.userRootPath,"/根目录"));
                                    gene_file=original_data.path;
                                    OriginalData.cansel_tree();
                                }else{
                                    appUtil.popbox ("您导入的数据不是fa文件，不符合Unigene数据的格式！！");
                                    return;
                                }

                            }else if("folder1" == type){
                                if(original_data.fileType=="folder"){
                                    $("#userName1").val(original_data.path.replace(jsonPath.userRootPath,"/根目录"));
                                    folder=original_data.path;
                                    OriginalData.cansel_tree();
                                }else{
                                    appUtil.popbox ("您导入的的数据不符合功能注释数据的格式,请重新导入！");
                                    return;
                                }
                            }
                        }
                        //$("#"+id).html("");
                    }
                }
            }
            OriginalData.cansel_tree =function(){
                $(".tree_show_id").hide();
            }
            //原始数据添加样品
            function _oriSelect(datas,type,filePath) {
                if(oneOrTwo=="single"){
                    Process.OriginalData.originalDataSin = Process.OriginalData.originalDataSin || {items: []};
                    Process.OriginalData.originalDataSin.size ? Process.OriginalData.originalDataSin.size : Process.OriginalData.originalDataSin.size = 0;
                    Process.OriginalData.originalDataSin.oldSize = Process.OriginalData.originalDataSin.size;
                }else if(oneOrTwo=="multiple"){
                    Process.OriginalData.originalDataMul = Process.OriginalData.originalDataMul || {items: []};
                    Process.OriginalData.originalDataMul.size ? Process.OriginalData.originalDataMul.size : Process.OriginalData.originalDataMul.size = 0;
                    Process.OriginalData.originalDataMul.oldSize = Process.OriginalData.originalDataMul.size;
                }
                var foldjQuery = $(els.fold);
                var filejQueryOne = $(els.file1);
                var filejQueryTwo = $(els.file2);
                var showData = $(".SampleInformation");
                var dataNew=[],num=0;

                switch (type) {//文件或者文件夹的选择类型
                    case "1"://选择文件夹
                        if(oneOrTwo=="single"){
                            if(JSON.stringify(datas)=="{}"){
                                appUtil.popbox ("您选择的文件夹不符合单端数据，请重新选择！");
                            }else{
                                for (var d in datas) {
                                    var filename = datas[d];
                                    if (filename) {
                                        //返回的数据可能出现文件颠倒的
                                        var name = filename.fileName1;
                                        if (name.indexOf("1.fq") < 0) {
                                            var name1 = filename.fileName1,
                                                path1 = filename.fq1,
                                                id1 = filename.fileid1,
                                                row1 = filename.raw_fq1;
                                            filename.fileName1 = name1;
                                            filename.fq1 = path1;
                                            filename.raw_fq1 = row1;
                                            filename.fileid1 = id1;
                                        }

                                        if (name != null) {
                                            if (!_isOriSelect(filename.fq1)) {
                                                /*datas[d].name = OriginalData._namePrefix(name);*/
                                                datas[d].name = OriginalData._namePrefix(name);
                                                dataNew[num]=datas[d];
                                                num++;
                                            }
                                        }
                                        /*if (!_isOriSelect(filename.fq1)) {
                                            var  num=parseInt(1000*Math.random());
                                            datas[d].name ="S0"+num;
                                        } else {
                                            datas.splice(d);
                                        }*/
                                    }
                                }
                                Process.OriginalData.originalDataSin.items.push.apply(Process.OriginalData.originalDataSin.items, dataNew);
                                $("#originalDataLi1").tmpl(Process.OriginalData.originalDataSin).appendTo(".SampleInformation");
                                OriginalData.cansel_tree();
                            }

                        }else if(oneOrTwo=="multiple"){
                            for (var d in datas) {
                                var filename = datas[d];
                                if (filename) {
                                    //返回的数据可能出现文件颠倒的
                                    var name = filename.fileName1;
                                    if (name.indexOf("1.fq") < 0) {
                                        var name1 = filename.fileName1,
                                            path1 = filename.fq1,
                                            id1 = filename.fileid1,
                                            row1 = filename.raw_fq1;
                                        filename.fileName1 = filename.fileName2;
                                        filename.fileName2 = name1;
                                        filename.fq1 = filename.fq2;
                                        filename.fq2 = path1;
                                        filename.raw_fq1 = filename.raw_fq2;
                                        filename.raw_fq2 = row1;
                                        filename.fileid1 = filename.fileid2;
                                        filename.fileid2 = id1;
                                    }
                                     if (name != null) {
                                         if (!_isOriSelect(filename.fq1)) {
                                             /*datas[d].name = OriginalData._namePrefix(name);*/
                                             datas[d].name = OriginalData._namePrefix(name);
                                             dataNew[num]=datas[d];
                                             num++;
                                         }/* else {
                                             datas.splice(d);
                                         }*/
                                     }
                                    /*var  num=parseInt(1000*Math.random());
                                    datas[d].name ="S0"+num;*/
                                }
                            }
                            Process.OriginalData.originalDataMul.items.push.apply(Process.OriginalData.originalDataMul.items, dataNew);
                            $("#originalDataLi").tmpl(Process.OriginalData.originalDataMul).appendTo(".SampleInformation");
                            OriginalData.cansel_tree();
                        }
                        break;
                    case "fastq1"://fastq1
                        //_get_file.apply(this, [filejQueryTwo, filejQueryOne, datas, showData, '2']);
                        _getFile(filejQueryOne, filejQueryTwo, datas, showData, type, filePath);
                        break;
                    case "fastq2"://fastq2
                        //_get_file.apply(this, [filejQueryOne, filejQueryTwo, treeObj, showData, '3']);
                        _getFile(filejQueryOne, filejQueryTwo, datas, showData, type, filePath);
                        break;
                    case "5":
                        var data = this.postJSON("/report/dge/getPath", {id: treeObj.getSelectedNodes()[0].id});
                        folder=treeObj.getSelectedNodes()[0].name;
                        break;
                    case "outdir":
                        Process.outputDir = treeObj.getSelectedNodes()[0].id;
                        Process.submitProcessStart();
                        removeFile();
                        break;
                    case "analysis-pid" :
                        $("#data-analysis-id").attr("data-pid", treeObj.getSelectedNodes()[0].id);
                        $("#data-analysis-id").val(treeObj.getSelectedNodes()[0].name);
                        break;
                    default:
                        break;
                }
            };
            //手动输入提交  添加联动效果,选择一个,另一个自动添加,并判断是否存在,不存在提示
            function _getFile(filejQueryOne, filejQueryTwo, datas, showData, type, filePath) {
                    if (datas.length == 0 && type == "fastq1") {
                        var filejQueryOther = filejQueryTwo.val().trim();
                        if (filejQueryOther == "") {
                            filejQueryOne.val("");
                            filejQueryOne.val(filePath);
                            appUtil.popbox ("没有找到对应的fq2文件");

                        } else {
                            var datas = [];
                            var file = {};
                            var file1 = filePath.split("/");
                            var file2 = filejQueryOther.split("/");
                            var fileName1 = file1[file1.length - 1];
                            var fileName2 = file2[file2.length - 1];
                            var gUID = appUtil.getUID();
                            file.fileName1 = fileName1;
                            file.fq1 = filePath;
                            file.raw_fq1 = filePath;
                            file.fileName2 = fileName2;
                            file.fq2 = filejQueryOther;
                            file.raw_fq2 = filejQueryOther;
                            file.ok = "false";
                            file.id = gUID;
                            datas[0] = file;
                            _getFile_from(datas);
                            filejQueryTwo.val("");
                        }

                    } else if (datas.length == 0 && type == "fastq2") {
                        var filejQueryOther = filejQueryOne.val().trim();
                        if (filejQueryOther == "") {
                            filejQueryTwo.val("");
                            filejQueryTwo.val(filePath);
                            appUtil.popbox ("没有找到对应的fq1文件");
                        } else {
                            var datas = [];
                            var file = {};
                            var file2 = filePath.split("/");
                            var file1 = filejQueryOther.split("/");
                            var fileName2 = file1[file1.length - 1];
                            var fileName1 = file2[file2.length - 1];
                            var gUID = appUtil.getUID();
                            file.fileName1 = fileName1;
                            file.fq1 = filejQueryOther;
                            file.raw_fq1 = filejQueryOther;
                            file.fileName2 = fileName2;
                            file.fq2 = filePath;
                            file.raw_fq2 = filePath;
                            file.ok = "false";
                            file.id = gUID;
                            datas[0] = file;
                            _getFile_from(datas);
                            filejQueryOne.val("");
                        }
                    } else {
                        _getFile_from(datas);
                    }
                }
            //提交选择的单个文件
            function _getFile_from(datas) {
                var items;
                if(oneOrTwo=="single"){
                    items = Process.OriginalData.originalDataSin.items;
                    for (var d in items) {
                        if ((datas[0].fileName1 == items[d].fileName1)) {
                            appUtil.popbox ("该fq文件已经存在！");
                            return;
                        }
                    }
                    datas[0].name =OriginalData._namePrefix(datas[0].fileName1);
                    Process.OriginalData.originalDataSin.items.push.apply(Process.OriginalData.originalDataSin.items, datas);
                    $(".SampleInformation").html("");
                    $("#originalDataLi1").tmpl(Process.OriginalData.originalDataSin).appendTo(".SampleInformation");
                }else if(oneOrTwo=="multiple"){
                    items = Process.OriginalData.originalDataMul.items;
                    for (var d in items) {
                        if ((datas[0].fq1 == items[d].fq1) && (datas[0].fq2 == items[d].fq2)) {
                            appUtil.popbox ("该fq文件已经存在！");
                            return;
                        }
                    }
                    datas[0].name =OriginalData._namePrefix(datas[0].fileName1);
                    Process.OriginalData.originalDataMul.items.push.apply(Process.OriginalData.originalDataMul.items, datas);
                    $(".SampleInformation").html("");
                    $("#originalDataLi").tmpl(Process.OriginalData.originalDataMul).appendTo(".SampleInformation");
                }
            }

            //用来关闭原始弹出框
            OriginalData.oriClosePop = function () {
                var nameStri = [];
                if($("#single").is(":checked")){
                    this.originalDataSin ? this.originalDataSin.size = this.originalDataSin.oldSize : null;
                }else if($("#multiple").is(":checked")){
                    this.originalDataMul ? this.originalDataMul.size = this.originalDataMul.oldSize : null;
                }

                $(".riCon").each(function () {
                    nameStri.push($(this).val());
                    for (var i in dataCfg) {
                        if (i == $(this).attr("content") && $(this).attr("content") != $(this).val()) {
                            var value = dataCfg[i];
                            dataCfg[$(this).val()] = value;
                            delete  dataCfg[i];
                        }
                    }
                });
                if (document.getElementById("a31c8e2c1feb4d509806c3a0bfbe2b9e")) {
                    this.removeOriginalDataList = [];
                    var oriData;
                    if($("#single").is(":checked")){
                        oriData=this.originalDataSin;
                    }else if($("#multiple").is(":checked")){
                        oriData=this.originalDataMul;
                    }
                    if (oriData && oriData.items) {
                        for (var i = 0; i < oriData.items.length; i++) {
                            if (!oriData.items[i].ok) {
                                oriData.items.splice(i, 1);
                                i--;
                            }
                            if (i == -1) {
                                continue;
                            }
                            oriData.items[i].name = nameStri[i];
                        }
                    }
                }
                //_itemsVSdataContent();
                $(".pop_up").remove();
            };
            return OriginalData;
        })();

        //差异分组选择模块代码
        var DegGroupingSelect = (function () {
            var DegGroupingSelect = new Object();
            _extends(DegGroupingSelect, Util);
            var els = {
                "popClick": "#degGroupingSelectPop",  //差异分组确定按钮id
                "popHtml": "#degGroupingSelectPopHtml", //差异分组弹出框html
                "firstOk": "#degSelectFirstOk", //差异分组弹出框中的第一个确定按钮id
                "con": ".differences input[name=Contrast_check]", //对照
                "exp": ".differences input[name=Experiment_check]", //实验
                "all": ".differences input", //左右input
                "group": "#E3A19A7D",  //组合div的id
                "groupLi": "#E3A19A7D div"//组合div的id
            };
            //差异分组弹出框弹出代码
            DegGroupingSelect.pop = function () {
                    if (!_isNullOfDataCfg()) {
                        DegGroupingSelect.degGoupSelect();
                    } else {
                        appUtil.popbox ("您还没有原始数据，无法进行其他操作");
                    }
            };
            //原始数据样品改变时
            DegGroupingSelect.onchangename = function (_this) {
                var oldName = $(_this).attr("content"),
                    newName = $(_this).val();
                if (detailCfg.Com) {
                    detailCfg.Com = _changeName(detailCfg.Com, oldName, newName);
                }
                if (detailCfg.Sep) {
                    detailCfg.Sep = detailCfg.Sep.replace(oldName, newName);
                }
            }
            /*DegGroupingSelect.onchangename = function (_this) {
                if($(_this).val()==""||$(_this).val()==null){
                    appUtil.popbox ("请输入样品名称!");
                }
                var oldName = $(_this).attr("content"),
                    newName = $(_this).val();
                //修改dataCfg中样本名
                for (var i in dataCfg) {
                    if (i == oldName && oldName != newName) {
                        var value = dataCfg[i];
                        dataCfg[newName] = value;
                        delete  dataCfg[i];
                    }
                }
                var oriData;
                if(oneOrTwo=="single"){
                    oriData=Process.OriginalData.originalDataSin;
                }else if(oneOrTwo=="multiple"){
                    oriData=Process.OriginalData.originalDataMul;
                }
                //修改差异分组中样本名
                for (var i = 0; i < oriData.items.length; i++) {
                    var item = oriData.items[i];
                    for (var name in item) {
                        if (name == 'name'&& item[name]==oldName) {
                            item[name]=newName;
                        }
                    }
                }
                if (detailCfg.Com) {
                    detailCfg.Com = _changeName(detailCfg.Com, oldName, newName);
                }
                if (detailCfg.Sep) {
                    detailCfg.Sep = detailCfg.Sep.replace(oldName, newName);
                }
            }*/

            function _changeName(data, oldName, newName) {
                var one = data.split("!"), dataOut = "";
                for (var i = 0; i < one.length; i++) {
                    dataOut += one[i].replace(oldName, newName) + (i == one.length ? "" : "!");
                }
                return dataOut;
            }
            //差异分组弹出框初始化
            DegGroupingSelect.degGoupSelect = function () {
                if (dataCfg.length == 0) {
                    throw "差异分组的配置文件没有初始化，请初始化在试？";
                    return;
                } else {
                    var gUID = this.getUID();
                    this.gUID = gUID;  //用来存放确定按钮的id
                    var jqueryTmplData = {"items": [], "gUID": gUID, "results": []};
                    var oriData;
                    if(oneOrTwo=="single"){
                        oriData=Process.OriginalData.originalDataSin.items;
                    }else if(oneOrTwo=="multiple"){
                        oriData=Process.OriginalData.originalDataMul.items;
                    }
                    for (var i = 0; i < oriData.length; i++) {
                        var item = oriData[i];
                        for (var name in item) {
                            if (name == 'name') {
                                var arr = jqueryTmplData.items;
                                arr[arr.length] = item[name];

                            }
                        }
                    }

                    _addJqueryTmplData(detailCfg.Com, arr, jqueryTmplData);
                    _addJqueryTmplData(detailCfg.Sep, arr, jqueryTmplData);
                    jqueryTmplData.items.sort(function (a, b) {
                        var a1 = parseInt(a.substring(1));
                        var b1 = parseInt(b.substring(1));
                        if (a1 < b1) {
                            return -1;
                        } else if (a1 > b1) {
                            return 1;
                        } else {
                            return 0;
                        }
                    });
                    $(".diffGroupModel").html("");
                    $("#degGroupingSelectPopHtml").tmpl(jqueryTmplData).prependTo(".diffGroupModel");//主流程第四步差异分组
                    this.selectOk();
                }
            };

            function _addJqueryTmplData(str, arr, jqueryTmplData) {
                if (str) {
                    var ssss1 = str.split("!");
                    for (var m = 0; m < ssss1.length; m++) {
                        if (_removeUserlessitem(arr, ssss1[m])) {
                            jqueryTmplData.results.push(ssss1[m]);
                        }
                    }
                }
            };

            //差异分组第一个确定按钮（添加按钮）
            DegGroupingSelect.selectOk = function () {
                $(els.firstOk).click(function () {
                    var contrasts = $(els.con);
                    var str1 = _addToArr(contrasts);
                    var experiments = $(els.exp);
                    var str2 = _addToArr(experiments);
                    var values = "";
                    if (str1.length == 0 || str2.length == 0) {
                        return;
                    }
                    if (str1.length == 1 && str2.length == 1) {
                        values = str1[0] + "," + str2[0];
                    } else {
                        for (var i = 0; i < str1.length; i++) {
                            (i == 0) ? (values += str1[i]) : (values += "," + str1[i]);
                        }
                        for (var i = 0; i < str2.length; i++) {
                            (i == 0) ? (values += ";" + str2[i]) : (values += "," + str2[i]);
                        }
                    }
                    var lis = $("#E3A19A7D div");
                    for (var i = 0; i < lis.length; i++) {
                        if ($(lis[i]).attr("data") == values) {
                            appUtil.popbox ("已有该差异分组，不能重复提交！");
                            return;
                        }
                    }
                    if($(els.group).children(".in").length == 0 ) {
                        $(els.group).append("<div data='"+values+"' class='text-center differ-group-div in'>" +
                        "<span>" + ((values.indexOf(";") > -1) ? (values.replace(/,/g, "_").replace(/;/g, "_vs_")) : (values.replace(/,/g, "_vs_")))+
                        "</span>&nbsp;&nbsp;<span><a href='javascript:void(0);' class='icon-close' onclick='DgeReport.PageStyle.callCode.removeLi(this);'><i class='fa fa-times'></i></a></span></div>");
                    }else{
                        $(els.group).find("div:first-child").before("<div data='"+values+"' class='text-center differ-group-div in'>" +
                        "<span>" + ((values.indexOf(";") > -1) ? (values.replace(/,/g, "_").replace(/;/g, "_vs_")) : (values.replace(/,/g, "_vs_")))+
                        "</span>&nbsp;&nbsp;<span><a href='javascript:void(0);' class='icon-close' onclick='DgeReport.PageStyle.callCode.removeLi(this);'><i class='fa fa-times'></i></a></span></div>");
                    }
                    $(els.all).attr("checked",false);
                    $(".differences").each(function () {
                        $(this).find("li input[disabled='disabled']").removeAttr("disabled");
                    });
                });
            };

            function _addToArr(oriArr) {
                var arr = [];
                for (var i = 0; i < oriArr.length; i++) {
                    if (oriArr[i].checked) {
                        arr[arr.length] = $(oriArr[i]).val();
                    }
                }
                return arr;
            };
            //差异分组第二个确定
            DegGroupingSelect.submit = function () {
                var lis = $(els.groupLi);
                detailCfg.Sep = "";//四个的分组
                detailCfg.Com = "";//两个的分组
                for (var i = 0; i < lis.length; i++) {
                    var datas = $(lis[i]).attr("data");
                    if (datas.indexOf(";") > -1) {
                        detailCfg.Sep ? (detailCfg.Sep += "!" + datas) : (detailCfg.Sep = datas);
                    } else {
                        detailCfg.Com ? (detailCfg.Com += "!" + datas) : (detailCfg.Com = datas);
                    }
                }
                if (!detailCfg.Com && !detailCfg.Sep) {
                    appUtil.popbox ("请您选择差异分组后提交任务。");
                    return false;
                }else{
                    return true;
                }
            };
            return DegGroupingSelect;
        })();

        //用来判断dataCfg中是否含有属性，也就是是否为空对象
        //true 表示 是空对象 false 表示非空对象
        function _isNullOfDataCfg() {
            if (dataCfg) {
                var i = 0;
                for (var key in dataCfg) {
                    i++;
                }
                if (i == 0) {
                    return true;
                } else {
                    return false;
                }
            } else {
                return true;
            }
        };

        //用来判断差异分组中的是否存在不存在的原始数据
        //返回 ture 表示不存在 也就是说差异分组中所有的样品都是合理的
        //返回false 表示存在  也就是说差异分组中存在某些样品是不合理的
        function _removeUserlessitem(arr, str) {
            var strarr = [];
            var ss = str.split(";");
            strarr.push.apply(strarr, ss[0].split(","));
            if (ss[1]) {
                strarr.push.apply(strarr, ss[1].split(","));
            }
            for (var i = 0, l = strarr.length; i < l; i++) {
                if (!_con(arr, strarr[i])) {
                    return false;
                }
            }
            return true;
        };
        function _con(arr, str) {
            for (var i = 0, l = arr.length; i < l; i++) {
                if (str == arr[i]) {
                    return true;
                }
            }
            return false;
        };

        //比较items 和dataContent的区别dataContent向items靠拢
        function _itemsVSdataContent() {
            var oriData;
            if(oneOrTwo==""){
                if(typeselect=="single"){
                    oriData=Process.OriginalData.originalDataSin;
                }else if(typeselect=="multiple"){
                    oriData=Process.OriginalData.originalDataMul;
                }
            }else if(oneOrTwo=="single"){
                oriData=Process.OriginalData.originalDataSin;
            }else if(oneOrTwo=="multiple"){
                oriData=Process.OriginalData.originalDataMul;
            }
            var items = oriData.items, parms = {};
            ;
            if (modeler(items) < modeler(dataContent)) {
                for (var i in dataContent) {
                    if (dataContent[i]) {
                        var count = 0;
                        for (var item in items) {
                            if (items[item].fileid1) {
                                parms[OriginalData._namePrefix(items[item].fileName1)] = items[item];
                            } else if (i == OriginalData._namePrefix(items[item].fileName1)) {
                                count++;
                            }
                        }
                        count == 0 ? delete dataContent[i] : "";
                    }
                }
            } else if (modeler(items) == modeler(dataContent)) {
                for (var item in items) {
                    if (item) {
                        var filename1 = OriginalData._namePrefix(items[item].fileName1);
                        if (items[item].fileid1) {
                            parms[filename1] = items[item];
                        }
                    }
                }
            } else if (modeler(items) > modeler(dataContent)) {
                for (var item in items) {
                    if (item) {
                        var count = 0, filename1 = OriginalData._namePrefix(items[item].fileName1);
                        if (items[item].fileid1) {
                            parms[filename1] = items[item];
                        } else {
                            for (var i in dataContent) {
                                if (i == filename1) {
                                    count++;
                                }
                            }
                            count == 0 ? parms[filename1] = items[item] : "";
                        }
                    }
                }
            }
            _getFilePath(parms);
        }

        //获得json的长度
        var modeler = function (obj) {

            var count = 0;
            if (obj && typeof obj === "object") {
                for (var ooo in obj) {
                    if (obj.hasOwnProperty(ooo)) {
                        count++;
                    }
                }
                return count;
            } else {
                return 0;
            }

        };
        //修改datacontent中的名字
        function _fromItemsToName(name) {
            var oriData;
            if(oneOrTwo==""){
                if(typeselect=="single"){
                    oriData=Process.OriginalData.originalDataSin;
                }else if(typeselect=="multiple"){
                    oriData=Process.OriginalData.originalDataMul;
                }
            }else if(oneOrTwo=="single"){
                oriData=Process.OriginalData.originalDataSin;
            }else if(oneOrTwo=="multiple"){
                oriData=Process.OriginalData.originalDataMul;
            }
            var json = oriData.items;
            for (var i in json) {
                var item = json[i];
                if (item.fileName1 == name) {
                    return item.name;
                }

            }
        }

        //样品文件相对路径生成配置文件dataContent.cfg  主流程提交按钮上面的dataContent里面的信息
        function _getFilePath(parms) {
            var json,
                num;
            if(oneOrTwo==""){
                if(typeselect=="single"){
                    num=1;
                    json = (parms == null ? Process.OriginalData.originalDataSin.items : parms), data = "";
                }else if(typeselect=="multiple"){
                    num=2;
                    json = (parms == null ? Process.OriginalData.originalDataMul.items : parms), data = "";
                }
            }else if(oneOrTwo=="single"){
                num=1;
                json = (parms == null ? Process.OriginalData.originalDataSin.items : parms), data = "";
            }else if(oneOrTwo=="multiple"){
                num=2;
                json = (parms == null ? Process.OriginalData.originalDataMul.items : parms), data = "";
            }

            for (var item in json) {
                if (item) {
                    if(num==1){
                        var file1 = json[item].fileName1;
                        if (_getidfromName(file1, json)) {
                            data += file1 + ":" + _getidfromName(file1, json)  + ",";
                        }
                    }else if(num==2){
                        var file1 = json[item].fileName1, file2 = json[item].fileName2;
                        if (_getidfromName(file1, json) && _getidfromName(file2, json)) {
                            data += file1 + ":" + _getidfromName(file1, json) + "," + file2 + ":" + _getidfromName(file2, json) + ",";
                        }
                    }

                }
            }
            if (data != null && data.length > 0) {
                var path = appUtil.postJSON("/report/dge/getFilePath", {data: data.substring(0, data.length - 1)});
                if (!jQuery.isEmptyObject(path)) {
                    for (var i in path.path) {
                        var one = path.path[i];
                        for (var o in one) {
                            if (o.indexOf("1.fq") > 0) {
                                var name = _fromItemsToName(o);
                                dataContent[name] = {
                                    'fq1': one[o] + '/' + o,
                                    'filename1': o  ,
                                    'raw_fq1': one[o] + '/' + o,
                                    'fq2': one[o] + '/' + (o.replace("1.fq", "2.fq")),
                                    'filename2': o.replace("1.fq", "2.fq"),
                                    'raw_fq2': one[o] + '/' + o.replace("1.fq", "2.fq")
                                };
                            }
                        }
                    }
                }
            }
        }

        function _getidfromName(name, json) {
            for (var i in json) {
                if (json[i].fileName1 == name) {
                    return json[i].fileid1;
                } else if (json[i].fileName2 == name) {
                    return json[i].fileid2;
                }
            }
        }

        //主流程点击提交按钮
        Process.submit = function (_this) {
            if (!appUtil.getPerm("5D6F471F-0022-4E13-B7D2-DB21C2DB540A")) {
                appUtil.popbox ("您没有使用的权限,请提升用户等级!");
                return;
            }
            var i = 0;
            for (var key in dataCfg) {
                i++;
            }
            if (i == 0) {
                appUtil.popbox ("原始数据为空，不可以提交任务。");
                return;
            }
            _getPageData();
            _getFilePath();
            if (!FDRdefined.trim() || !FCdefined.trim()) {
                appUtil.popbox ("FDR或FC为空，不可以提交任务。");
                return;
            }
            if (!detailCfg.Com && !detailCfg.Sep) {
                appUtil.popbox ("请您选择差异分组后提交任务。");
                return;
            }

            this.startStep = "S1";
            var isSubmitProcess = false;
            var qphredBool="33";
            $(".processCustom").find("input[name='Qphred']").each(function(){
                if($(this).is(":checked")==true){
                    qphredBool = $(this).val();
                }
            });
            if (this.isBlank) {
                if (!(dataCfg && (SepDefined || ComDefined))) {
                    isSubmitProcess = false;
                } else {
                    isSubmitProcess = true;
                }
            } else {
                var ref = new RegExp('^[A-Za-z]+$', 'g');
                if(!ref.test($("#titInput").val())){
                    appUtil.popbox ("您输入的'物种名称'信息不正确，请您调整后再提交任务。");
                    return ;
                }
                //校验FDR和FC的值
                var FDRValue=/*/^\d+(?=\.{0,1}\d+$|$)/;00000*//^0\.[0-9]{1,5}?$/;/*/^\d+(?=\.{0,1}\d+$|$)/;*/
                if(document.getElementById("FDRVal")){
                    if(! FDRValue.test($("#FDRVal").val())){
                        appUtil.popbox ("您输入的FDR的值不正确,请您调整后再提交任务。");
                        return ;
                    }
                }
                if(document.getElementById("FCVal")){
                    var FCValue=  new RegExp('^[0-9]+$', 'g');
                    if( !FCValue.test($("#FCVal").val())){
                        appUtil.popbox ("您输入的差异筛选倍数阈值不正确,请您调整后再提交任务。");
                        return ;
                    }
                }
                var dataCfgStr = JSON.stringify(dataCfg);
                var dataCfgNextStr = JSON.stringify(this.dataCfgNext);
                if (detailCfg.Sep != SepDefined ||
                    detailCfg.Com != ComDefined ||
                    detailCfg.fold != FCdefined ||
                    detailCfg.FDR != FDRdefined ||
                    annoDatabseTypeNew != this.annoDatabseType ||
                    dataCfgStr != dataCfgNextStr ||
                    $("#titInput").val() != $("#speciesName").val() ||
                    _getRefSpecies()||
                    qphredBool != "33" ||
                    $("#7C57F3B43FC04B6A9CE3DD91FDCDB534").val() != "fr-unstranded" ||
                    $("#7C8CEE4D2E264AD2B3D67601B0DA1A5A").val()!="NONE"
                ) {
                    _isStepStart.apply(this, null);
                    isSubmitProcess = true;//可以提交流程定制任务
                } else {
                    isSubmitProcess = false;
                    appUtil.popbox ("您输入的信息不完整或者功能模块参数无变化，请您调整后再提交任务。");
                    return;
                }
            }
            if (isSubmitProcess) {
                var data = new Date(),
                    dataCfgStr = JSON.stringify(dataCfg),
                    dataCfgNextStr = JSON.stringify(this.dataCfgNext),
                    outPathDate = "dge_" +
                        $("#titInput").val() + "_" + data.getFullYear() + ((data.getMonth() + 1) < 10 ? ("0" + (data.getMonth() + 1)) : (data.getMonth() + 1)) +
                        (data.getUTCDate() < 10 ? "0" + data.getUTCDate() : data.getUTCDate()) +
                        (data.getHours() < 10 ? "0" + data.getHours() : data.getHours()) +
                        (data.getMinutes() < 10 ? "0" + data.getMinutes() : data.getMinutes()) +
                        (data.getSeconds() < 10 ? "0" + data.getSeconds() : data.getSeconds());
                if (dataCfgStr != dataCfgNextStr || $("#titInput").val() != $("#speciesName").val() || qphredBool != "33" || $("input[name=Lib_type]:checked").val() != "fr-unstranded" ) {
                    //新项目
                    projectCode = outPathDate;
                    _typeproject = "1.0";
                } else {
                    //新版本
                    var josnType = appUtil.postJSON("/report/dge/typeproject", {projectCode: $("#projectCode").val() == "" ? outPathDate : $("#projectCode").val()});
                    _typeproject = josnType.version == null ? "1.0" : josnType.version;
                    projectCode = $("#projectCode").val();
                }
                var jqueryTmplData = {}, data = new Date(), dataString = "",
                    com = "", sep = "", LinuxPath,
                    wuzhong = $("#annoDatabseType").val(),
                    LinuxPath = _projectPath = "/biomarker_project/" + projectCode + "/v" + _typeproject + "/" + outPathDate + "/";
                if (detailCfg.Sep != null && detailCfg.Sep != "") {
                    var string = detailCfg.Sep.split("!");
                    for (var i = 0; i < string.length; i++) {
                        var temp = string[i].split(";");
                        for(var ii in temp){
                            sep += temp[ii].replace(/,/g, "_") ;
                            ii == temp.length-1 ? "" : sep += "_vs_";
                        }
                        i == string.length -1 ? "" : sep+= ",";
                    }
                }
                if (detailCfg.Com != null && detailCfg.Com != "") {
                    var string = detailCfg.Com.split("!");
                    for (var i = 0; i < string.length; i++) {
                        com += string[i].replace(",", "_vs_") + ",";
                    }
                }
                var jsonPath=appUtil.postJSON("/report/util/getUserPath");
                for (var i in dataCfg) {
                    for (var item in dataCfg[i]) {
                        if (item == 'fq1' || item == 'fq2') {
                            if((dataCfg[i][item]+"").contains(jsonPath.userRootPath)){
                                dataString += "<p>" + i + "\t" + item + "\t" + dataCfg[i][item].replace(jsonPath.userRootPath,"/根目录") + "</p>";
                            }else{
                                var str="",filename=dataCfg[i][item],length=filename.split("/").length;
                                str+="/根目录/"+filename.split("/")[length-2]+"/"+filename.split("/")[length-1];
                                dataString += "<p>" + i + "\t" + item + "\t" + str + "</p>";
                            }
                        }
                    }
                }

                jqueryTmplData.shuchu = "/根目录" + LinuxPath;
                jqueryTmplData.beishu = FCdefined;
                jqueryTmplData.FDR = FDRdefined;
                jqueryTmplData.sepcom = (sep == "" ? "" : (sep + ",")) + com.substring(0, com.length - 1);
                var  versions=$("#11523B85157341A4A31F118243F568B2").find("input[type='radio']");//组装版本
                if(versions.is(":checked")){
                    jqueryTmplData.refSpVer=$(versions).val();
                }
                jqueryTmplData.refSp =$("#7C8CEE4D2E264AD2B3D67601B0DA1A5A").find(".checked").find("a").text();//参考物种
                jqueryTmplData.projectName = $("#project_ID").val().trim();
                jqueryTmplData.folder = folder;
                jqueryTmplData.gene_file = gene_file;
                //miaorong修改
                var COMString = "",
                    refSp = $("#7C8CEE4D2E264AD2B3D67601B0DA1A5A").find(".checked").find("a").text(),
                    refSpData = jqueryTmplData.refSpVer,qphredBool ="33";
                if(refSp != "NONE" ){
                    for(var items in Process.referenceSpecies){
                        if(items.indexOf(refSp) >= 0){
                            if(Process.referenceSpecies[items].split("\t")[2] == refSpData){
                                var detail = Process.referenceSpecies[items].split("\t");
                                detailCfg.Preparation = detail[4].split(":")[1];
                                detailCfg.Genome_unigene = detail[5].split(":")[1];
                                detailCfg.GFF = detail[6].split(":")[1];
                            }
                        }
                    }
                }
                var oriData;
                if(oneOrTwo==""){
                    if(typeselect=="single"){
                        oriData=Process.OriginalData.originalDataSin;
                    }else if(typeselect=="multiple"){
                        oriData=Process.OriginalData.originalDataMul;
                    }
                }else if(oneOrTwo=="single"){
                    oriData=Process.OriginalData.originalDataSin;
                }else if(oneOrTwo=="multiple"){
                    oriData=Process.OriginalData.originalDataMul;
                }
                for (var i = 0; i < oriData.items.length; i++) {
                    var items = oriData.items[i];
                    for (var item in items) {
                        if (item == 'name') {
                            COMString += items[item] + ",";
                        }
                    }
                }
                if (detailCfg.Com == "" || detailCfg.Com == null) {
                    delete  detailCfg.Com;
                }
                if (detailCfg.Sep == "" || detailCfg.Sep == null) {
                    delete detailCfg.Sep;
                }
                $(".processCustom .panel-body").find("input[name='Qphred']").each(function(){
                    if($(this).is(":checked")){
                        qphredBool = $(this).val();
                    }
                });
                var wuzhong = $("#annoDatabseType").val();
                for(var item in annotation){
                    if(item == wuzhong ){
                        detailCfg.Cog = annotation[item].Cog;
                        detailCfg.nr = annotation[item].nr;
                        detailCfg.Kegg = annotation[item].Kegg;
                        detailCfg.Swissprot = annotation[item].Swissprot;
                    }
                }
                dataCfg.Qphred = qphredBool;
                detailCfg.fold = FCdefined;
                detailCfg.Project_name = $("#project_ID").val();
                detailCfg.Customer_info = "百迈克生物科技有限公司";
                detailCfg.FDR = FDRdefined;
                detailCfg.Project_id = projectCode;
                detailCfg.Project_key = $("#titInput").val();

                if(parameter=="wu"){
                    detailCfg.Preparation = folder;
                    detailCfg.Genome_unigene = gene_file;
                    delete detailCfg.GFF;
                }

                console.log("提交的detail配置文件：" + JSON.stringify(detailCfg));
                console.log("提交的data配置文件：" + JSON.stringify(dataCfg));
                $('.submitTextDataModel').html("");
                $("#submitTextDataHtml").tmpl(jqueryTmplData).prependTo('.submitTextDataModel');
                $("#originalDatainput").html(dataString);
                if(parameter=="you"){
                    $(".wuzhong").show();
                    $(".dataImport").hide();
                }else if(parameter=="wu"){
                    $(".wuzhong").hide();
                    $(".dataImport").show();
                }
                //miaorong修改
            }
        };
        //比较时候修改了原始数据框
        function _comparisonData(dataCfg, originalDataLis) {
            var key = [], list = "";
            for (var item in dataCfg) {
                if (item) {
                    key.push(item);
                }
            }
            for (var item in originalDataLis) {
                if (item == 'name') {
                    if (!key.indexOf(originalDataLis[item])) {
                        return false;
                    }
                }
            }
            if (dataCfg.length == originalDataLis.length) {
                return true;
            } else {
                return false;
            }
        }
        //比较参考物种
        function _getRefSpecies (){
            var refSp,refSpVer;
            var  versions=$("#11523B85157341A4A31F118243F568B2").find("input[type='radio']");//组装版本
            if(versions.is(":checked")){
                refSpVer=$(versions).val();
            }
            refSp =$("#7C8CEE4D2E264AD2B3D67601B0DA1A5A").find(".checked").find("a").text();//参考物种
            if(refSp == "NONE"){
                //没有改变
                return false;
            }else{
                if(detailCfg.Preparation.indexOf(Process.refSpeciesData )>0 && detailCfg.Genome_unigene.indexOf(Process.refSpeciesData)>0
                    && detailCfg.GFF.indexOf(Process.refSpeciesData)>0 ){
                    return false;
                }else{
                    return true;
                }
            }
        }
        //获得差异分析的所有值

        function _getPageData() {
            /*FDRdefined = $("*[name=FDR]").val();
             FCdefined = $("*[name=FC]").val();*/

            //FDR自定义数据读取
            if($("input[name=FDR]:checked").val()=="Customer"){
                //$("#inputFDR").show();
                FDRdefined = $("input[name=FDR]").val();
            }else{
                //$("#inputFDR").hide();
                FDRdefined = $("input[name=FDR]:checked").val();
            }
            if($("input[name=FC]:checked").val()=="Customer"){
                //$("#inputFC").show();
                FCdefined=$("input[name=FC]").val();
            }else{
                //$("#inputFDR").hide();
                FCdefined=$("input[name=FC]:checked").val();
            }

            if (FDRdefined <= 0 || FDRdefined > 0.1) {
                appUtil.popbox ("请输入正确的FDR数值(0, 0.1]");
                throw "";
            }
            if (FCdefined < 0 || FCdefined > 64) {
                appUtil.popbox ("请输入正确的FC数值[0, 64]");
                throw "";
            }
            if (dataCfg) {
                var asdf = 1;
                for (var o in dataCfg) {
                    asdf++;
                }
                if (asdf < 2) {
                    para_K_covDefined = "1";
                } else if (2 <= asdf < 16) {
                    para_K_covDefined = "2";
                } else {
                    para_K_covDefined = "3";
                }
            }
            annoDatabseTypeNew = $("#annoDatabseType").val();
        };

        //判断流程定制式从哪一步开始运行
        function _isStepStart() {
            var dataCfgStr = JSON.stringify(dataCfg),
                dataCfgNextStr = JSON.stringify(this.dataCfgNext),
                titInput = $("#titInput").val(),
                speciesName = $("#speciesName").val(),qphredBool="33";
            $(".processCustom .panel-body").find("input[name='Qphred']").each(function(){
                if($(this).is(":checked")){
                    qphredBool = $(this).val();
                }
            });
            if(dataCfgStr != dataCfgNextStr||
                $("#titInput").val() != $("#speciesName").val() ||
                !_getRefSpecies() ||
                qphredBool != "33" ||
                $("input[name=Lib_type]:checked").val() != "fr-unstranded" ||
                annoDatabseTypeNew != $("#annoDatabseType").val()
            ){
                this.startStep = "S1";
            }else if (detailCfg.Sep != SepDefined ||
                detailCfg.Com != ComDefined ||
                detailCfg.fold != FCdefined ||
                detailCfg.FDR != FDRdefined ) {
                this.startStep = "4";
            }
        };
        //开始提交流程定制任务任务
        Process.submitProcessStart = function () {
            var proName=$("#project_ID").val().trim();//项目名称
            var patrn = new RegExp("^[\u4e00-\u9fa5a-zA-Z0-9_]+$");
            if(proName==""){
                appUtil.popbox ("项目名称不允许为空，请输入项目名称！");
                return;
            }
            if (!patrn.test(proName)) {
                appUtil.popbox ("项目名称只能由中文、字母和数字组成！");
                return;
            }
            var COMString = "",
                num,
                oriData,refSpData,
                refSp = $("#7C8CEE4D2E264AD2B3D67601B0DA1A5A").find(".checked").find("a").text(),qphredBool ="33";
                var  versions=$("#11523B85157341A4A31F118243F568B2").find("input[type='radio']");//组装版本
                if(versions.is(":checked")){
                    refSpData=$(versions).val();
                }
            if(refSp != "NONE" ){
                for(var items in Process.referenceSpecies){
                    if(items.indexOf(refSp) >= 0){
                        if(Process.referenceSpecies[items].split("\t")[2] == refSpData){
                            var detail = Process.referenceSpecies[items].split("\t");
                            if(parameter=="you"){
                                detailCfg.Preparation = detail[4].split(":")[1];
                                detailCfg.Genome_unigene = detail[5].split(":")[1];
                                detailCfg.GFF = detail[6].split(":")[1];
                            }
                        }
                    }
                }
            }
            if(oneOrTwo==""){
                if(typeselect=="single"){
                    num=1;
                }else if(typeselect=="multiple"){
                    num=2;
                }
            }else if(oneOrTwo=="single"){
                num=1;
            }else if(oneOrTwo=="multiple"){
                num=2;
            }
            if(num==1){
                oriData=Process.OriginalData.originalDataSin;
            }else if(num==2){
                oriData=Process.OriginalData.originalDataMul;
            }
            for (var i = 0; i < oriData.items.length; i++) {
                var items = oriData.items[i];
                for (var item in items) {
                    if (item == 'name') {
                        COMString += items[item] + ",";
                    }
                }
            }

            if (detailCfg.Com == "" || detailCfg.Com == null) {
                delete  detailCfg.Com;
            }
            if (detailCfg.Sep == "" || detailCfg.Sep == null) {
                delete detailCfg.Sep;
            }
            $(".processCustom .panel-body").find("input[name='Qphred']").each(function(){
                if($(this).is(":checked")){
                    qphredBool = $(this).val();
                }
            });
            var wuzhong = $("#annoDatabseType").val();
            for(var item in annotation){
                if(item == wuzhong ){
                    detailCfg.Cog = annotation[item].Cog;
                    detailCfg.nr = annotation[item].nr;
                    detailCfg.Kegg = annotation[item].Kegg;
                    detailCfg.Swissprot = annotation[item].Swissprot;
                }
            }
            dataCfg.Qphred = qphredBool;
            detailCfg.fold = FCdefined;
            detailCfg.Project_name = proName;
            detailCfg.Customer_info = "百迈克生物科技有限公司";
            detailCfg.FDR = FDRdefined;
            detailCfg.Project_id = projectCode;
            detailCfg.Project_key = $("#titInput").val();
            if(parameter=="wu"){
                detailCfg.Genome_unigene=gene_file;
                detailCfg.Preparation=folder;
                delete detailCfg.GFF;
            }
            console.log("提交的detail配置文件：" + JSON.stringify(detailCfg));
            console.log("提交的data配置文件：" + JSON.stringify(dataCfg));
            var json={data:dataCfg,detail:detailCfg};
            /*this.postJSONAsync("/report/util/cfgFileCom",{json:JSON.stringify(json),projectType:"dge"},function(msg){
                var message=msg.power;
                fileSize=msg.fileSize;
                fileSize=parseFloat(fileSize)/2;
                var userLevel=msg.UserLevel[0].roleEN;
                console.info(message);
            });*/

            this.postJSONAsync("/report/dge/submitProcess",
                {
                    detail: JSON.stringify(detailCfg),
                    data: JSON.stringify(dataCfg),
                    outPath: _projectPath,
                    step: this.startStep,
                    projectCode: projectCode,
                    speciesName: $("#titInput").val(),
                    projectId: $("#projectId").val(),
                    type: _typeproject,
                    projectName: "表达谱分析项目",
                    projectIdentification: ($("#projectIdentification").val() == null ? "" : $("#projectIdentification").val()),
                    fileSize:fileSize
                }, _subCallback1());
            Dge.PageStyle.callCode.closePop();
        };
        //提交任务后的回调函数
        function _subCallback1() {
            return function (data) {
                if (data.error != "") {
                    appUtil.popbox (data.error);
                }else if (data.jobName != "") {
                    if (confirm("您的任务已经提交，任务名为 " + data.jobName + "版本为 v" + data.version + "，可以在我的任务进行查看。\r\n主流程运行时间较长，您可以查看云平台其他功能，关闭窗口不影响任务进行。\r\n         确定需要跳转到<我的项目>页面么?")) {
                        window.open("/project/project/showmyproject?permDataId=2F2FC7F7-537E-498F-B53F-5D5AE6D722E7");
                    }
                } else {
                    Process.warnTipModal("任务提交失败，请重试！");
                }
            };
        };
        /*原始数据导入校验*/
        function judgeSamplesName(){
            var patrn = new RegExp("^[a-zA-Z][0-9a-zA-Z_]{0,}$");//样品名称只能由数字、字母和下划线(\"_\")组成，且必须以字母开头
            var samplesDiv=$("#a31c8e2c1feb4d509806c3a0bfbe2b9e");
            if(samplesDiv.find("li").length==0){
                appUtil.popbox ("您还没有导入原始数据，请导入原始数据！");
                return;
            }
            var sampleDiv=samplesDiv.find("div.col-sm-6");
            for (var i=0;i<sampleDiv.length;i++){
                var sampleName=$(sampleDiv[i]).find("input.samples-input").val().trim();
                var sampleOther=$(sampleDiv[i]).siblings("div.col-sm-6");
                if(!(patrn.test(sampleName))){
                    appUtil.popbox ("样品名称只能由数字、字母和下划线(\"_\")组成，且必须以字母开头！");
                    return ;
                }
                if(sampleName.length>7){
                    appUtil.popbox ("样品名称不能大于7个字符");
                    return ;
                }
                for(var j=0;j<sampleOther.length;j++){
                    var sampleOtherName=$(sampleOther[j]).find("input.samples-input").val().trim()
                    if(sampleName==sampleOtherName){
                        appUtil.popbox ("样品名称不能重复");
                        return ;
                    }
                }
            }
            return true;
        }


        Process.nextSteps=function (){
            $("#formProcess").steps({
                bodyTag: "fieldset",
                onStepChanging: function (event, currentIndex, newIndex) {
                    switch(currentIndex){
                        //选择有参还是无参
                        case 0:
                            $(".clearfix a[href='#finish']").hide();
                            var a=Process.OriginalData.projectInput();
                            var b=Process.OriginalData.focusInput();
                            if(a&&b){
                                return true;
                            }else{
                                appUtil.popbox ("物种名称或项目名称输入有误！");
                                return false;
                            }
                            break;
                        //原始数据导入
                        case 1:
                            $(".clearfix a[href='#finish']").hide();
                            if(judgeSamplesName()){
                                Process.OriginalData.submit();
                                Process.annoSpeciesNext();//初始化差异分组
                            }else{
                                return;
                            }
                            break;
                        //参考物种
                        case 2:
                            $(".clearfix a[href='#finish']").hide();
                            if(Process.submitReferSpecies()){
                                //DegGroupingSelect.pop();
                                return true;
                            }else{
                                return false;
                            }
                            break;
                        //差异分组选择
                        case 3:
                            $(".clearfix a[href='#finish']").hide();//隐藏finish
                            if(Process.differentialNext()){
                                Process.submit();
                                $(".clearfix a[href='#finish']").show();//隐藏finish
                                return true;
                            }else{
                                return false;
                            }
                            break;
                        //参数确认
                        case 4:
                            /*$(".clearfix a[href='#finish']").show();
                             if(judgeSamplesName()){
                             OriginalData.submit();
                             Process.annoSpeciesNext();//初始化差异分组
                             }else{
                             return;
                             }*/
                            break;
                    }
                    if (currentIndex > newIndex) {
                        return true;
                    }

                    // Forbid suppressing "Warning" step if the user is to young
                    if (newIndex === 3 && Number($("#age").val()) < 18) {
                        return false;
                    }

                    var form = $(this);

                    // Clean up if user went backward before
                    if (currentIndex < newIndex) {
                        // To remove error styles
                        $(".body:eq(" + newIndex + ") label.error", form).remove();
                        $(".body:eq(" + newIndex + ") .error", form).removeClass("error");
                    }

                    // Disable validation on fields that are disabled or hidden.
                    form.validate().settings.ignore = ":disabled,:hidden";

                    // Start validation; Prevent going forward if false
                    return form.valid();
                },
                onStepChanged: function (event, currentIndex, priorIndex) {
                    // Suppress (skip) "Warning" step if the user is old enough.
                    if (currentIndex === 2 && Number($("#age").val()) >= 18) {
                        $(this).steps("next");
                    }

                    // Suppress (skip) "Warning" step if the user is old enough and wants to the previous step.
                    if(currentIndex==1){
                    }
                    if(currentIndex==0||currentIndex==1||currentIndex==2||currentIndex==3){
                        $(".clearfix a[href='#finish']").hide();
                    }
                },
                onFinishing: function (event, currentIndex) {

                    Process.submitProcessStart();
                },
                onFinished: function (event, currentIndex) {
                    var form = $(this);

                    // Submit form input
                    form.submit();
                }
            }).validate({
                errorPlacement: function (error, element) {
                    element.before(error);
                },
                rules: {
                    confirm: {
                        equalTo: "#password"
                    }
                }
            });
        }

        Process.OriginalData = OriginalData;
        Process.DegGroupingSelect = DegGroupingSelect;
        return Process;
    })();
    //数据挖掘模块  页面上拆分成3部分,js代码不变,数据内容添加uuid
    DataMining = (function () {
        var DataMining = new Object();
        DataMining.sequenceSamples = "";
        _extends(DataMining, Util);

        var els = {
            "degSelect": "#degGroupingInit", //差异分组下拉框
            "annoSearch": "#annoSearchKeyword", //注释搜索按钮id
            "keywords": "#annoKeywords", //注释搜索关键字
            "searchType": "#js-search-keyword-type", //搜索类型
            "tabUl": ".content_tab_title", //数据挖掘标签ul id
            "tabNumber": ".content_tab_title li", //数据挖掘标签头部
            "MiningHtml": "#dataMiningContent", //数据挖掘的html模板
            "MiningTrendChart": "#dataMiningTrendChart",//基因表达趋势图
            "MiningCon": ".rna_no_ref_report_tab", //存放模板的div id
            "table": "#tableDataShow", // 用来显示数据挖掘中的表格数据
            "degResult": "#degGroupGetDataByType", //差异结果查看的提交按钮id
            "degResultSelect": "select[name=degGrouping]", //获得选择的差异结果
            "degResultRegulated": "input[name=degGroupingRegulated]", //用来获取差异结果查看的上下条的类别
            "seqAna": ".js-dataMining-sequence-analysis", //序列分析
            "seqDown": ".js-dataMining-sequence-download", //序列下载
            "listAna": ".js-dataMining-list-download", //列表下载
            "sortAnalysis": ".js-dataMining-sort-analysis",//排序筛选
            "analysis": "#saveFileToLook", //文件分析查看按钮
            "name": "#data-analysis-name",
            "pid": "#data-analysis-id",
            "proteinInteraction": ".js-dataMining-protein-interaction" //蛋白互作图提交
        };

        //序列搜索添加条件选择--反向互补搜索
        /*$(els.searchType).change(function () {
            if ($(this).val() == "SequenceSearch") {
                $("#idSequenceSearch").html("<input type='checkbox' style='margin: -1px 3px 0px 0px;' value='反向互补搜索' id='sequenceCheckBox'>反向互补搜索");
            } else {
                $("#idSequenceSearch").text("从所有样品注释表达量总表中检索：");
            }
        });*/
        //序列搜索添加条件选择--反向互补搜索
        $(els.searchType).change(function () {
            /*if (!appUtil.getPerm("E6748DE9-0EAC-4F37-96B0-3BE51745C373")) {
                appUtil.popbox ("您没有使用的权限,请提升用户等级!");
                return;
            }*/
            if ($(this).val() == "SequenceSearch") {
                $("#idSequenceSearch").hide();
                $("#idSequenceSearch").parent(".form-group").find(".checkbox").show();

            } else {
                $("#idSequenceSearch").show();
                $("#idSequenceSearch").parent(".form-group").find(".checkbox").hide();
            }
        });

        //显示数据挖掘内容
        function _showContentStyle(idHtml) {
            if (idHtml) {
                $("#" + idHtml).show();
            } else {
                $("#dawjdoaijfisndvihzvlarng").show();
            }
            //表达基因趋势图显示内容
            $(".problem_area").each(function (b) {
                var id = this.id;
                $(this).mouseover(function (d) {
                    $("." + id).siblings(".problem_content").hide();
                    $("." + id).show();
                    $("." + id).css({
                        weight: "14px",
                        position: "fixed",
                        left: d.pageX + "px",
                        top: d.pageY + "px",
                        opacity: "0.8"
                    }).show();
                }).mouseout(function () {
                    $("." + id).hover(function () {
                        $("." + id).siblings(".problem_content").hide();
                        $("." + id).show();
                    }, function () {
                        $("." + id).hide();
                    });
                    $("." + id).hide();
                })
            })
        }

        //初始化数据挖掘
        DataMining.init = function () {
            if (Process.detailCfgNext) {
                var seps="",Coms="",arr = [], str="", i, l, degGroupSelectArr = [];
                if(("Sep" in Process.detailCfgNext)){
                    seps = Process.detailCfgNext.Sep;
                }
                if(("Com" in Process.detailCfgNext)){
                    for(var item in Process.detailCfgNext){
                        if(item == "Com"){
                            Coms += Process.detailCfgNext[item];
                        }
                    }
                }
                if (seps) {
                    arr.push.apply(arr, seps.split("!"));
                }
                if (Coms) {
                    arr.push.apply(arr, Coms.split("!"));
                }
                if (!arr.length) {
                    $(els.degSelect).append("<option value=\"NONE\">NONE</option>");
                }
                for (i = 0, l = arr.length; i < l; i++) {
                    (arr[i].indexOf(";") > -1) ? (str = arr[i].replace(/,/g, "_").replace(/;/g, "_vs_")) : (str = arr[i].replace(/,/g, "_vs_"));
                    degGroupSelectArr.push(str);
                    $(els.degSelect).append("<option value=\"" + str + "\">" + str + "</option>");
                }
                this.degGroupSelectArr = degGroupSelectArr;
                $("ul.FPKMValue li").find("input[name='hierarchicalClustering']").attr("disabled", "disabled");
            } else {
                throw "数据初始化出错,请重试！";
            }
        };

        //给数据挖掘页面添加tab标签
        var taskCount = 1;
        //数据挖掘内容部分显示方式,type:作图结果显示方式, num: 表达趋势图方式true,filp:蛋白互作图,button:排序筛选是否可选 ; page :翻页设置
        function _addTabToElement(data, tabName, type, num, filp, button, table, page) {
            if (countTab(table)) {
                return;
            }
            tabName = tabName + taskCount++;
            var gUID = appUtil.getUID();
            var asd = {"gUID": gUID, "type": type};
            if (type == 2) {
                var value = $(els.degSelect).val();
                asd.deg = value;
            }
            asd.filp = filp == true ? true : false;//表达venn图沉默蛋白胡作图
            asd.page = page == true ? true : false;//翻页设置
            if (num) {
                $("#" + table + " " + els.tabUl).append("<li class=\"active\" onclick=\"DgeReport.DataMining.PageView.Switch(this);\" Relation=\"" + gUID + "\"><a >" + tabName + "&nbsp;<span onclick=\"DgeReport.DataMining.PageView.removeThisTab(this,'" + table + "');\">X</span></a></li>");

                $(els.MiningTrendChart).tmpl(asd).appendTo($("#" + table + " " + els.MiningCon));
            } else {
                $("#" + table + " " + els.tabUl).append("<li class=\"active\" onclick=\"DgeReport.DataMining.PageView.Switch(this);\" Relation=\"" + gUID + "\"><a href=\"#" + gUID + "\" >" + tabName + "&nbsp;<span onclick=\"DgeReport.DataMining.PageView.removeThisTab(this,'" + table + "');\">X</span></a></li>");

                $(els.MiningHtml).tmpl(asd).appendTo($("#" + table + " " + els.MiningCon));
            }
            //排序筛选沉默
            if (button == null || !button) {
                $("#" + table + " .searchResult[relation='" + gUID + "']").find(".js-dataMining-sort-analysis").attr("style", "display:none");
                $("#" + table + " .searchResult[relation='" + gUID + "']").find(".js-dataMining-sort-analysis").attr("class", "");
            }
            PageView.Switch($("#" + table + " li[Relation=" + gUID + "]").get(0)); //添加按钮切换事件
            PageStyle.callCode.zuotu(); //添加 做图转换事件
            (function () {
                $(".search_rna a").hover(function () {
                    $(this).siblings(".problem_content").show();
                }, function () {
                    $(this).siblings(".problem_content").hide();
                });
                $(".problemtit a").hover(function () {
                    $(this).siblings(".problem_content").show();
                }, function () {
                    $(this).siblings(".problem_content").hide();
                });
                /*DataMining.seqDownload(table);
                DataMining.listDownload(table);
                DataMining.listAna(table);*/
            })();
            return gUID;
        };

        //验证该页面能存放的最大标签,默认值为6个，可以手动修改该默认值
        function countTab(table) {
            var length = $("#" + table + " " + els.tabNumber).length;
            if (length >= _maxTabCount) {
                appUtil.popbox ("对不起，该页面只能显示" + _maxTabCount + "个标签，请您删除无用标签后在进行操作。");
                return true;
            }
            return false;
        };
        DataMining.GeneMiningCheck=function(_this){
            var str=$("#js-search-keyword-type").val();
            var key=$("#annoKeywords").val();
            if(key==""){
                appUtil.popbox ("搜索信息不能为空，请重新输入!");
                return false;
            }else{
                var reg;
                if(str=="annoKey"){//功能
                    reg=new RegExp("^[A-Za-z0-9]+$");
                    if(!reg.test(key)){
                        appUtil.popbox ("功能注释只能有字母和数字，请重新输入!");
                        return false;
                    } else {
                        if(key.length>20||key.length<0) {
                            appUtil.popbox ("功能注释只能输入20个字符，请重新输入!");
                            return false;
                        }
                        return true;
                    }
                }
                if(str=="idKey") {//id
                    /*reg = new RegExp("^[A-Za-z0-9]+$");
                    if (!reg.test(key)) {
                        appUtil.popbox ("基因ID只能有字母和数字，请重新输入!");
                        return false;
                    } else {
                        return true;
                    }*/
                    return true;
                }
                if(str=="SequenceSearch"){//序列片段
                    reg = new RegExp("^[TACGtacg]+$", "g");
                    if(!reg.test(key)){
                        appUtil.popbox ("序列片段只能有大写字符TACG组成，请重新输入!");
                        return false;
                    } else {
                        if(key.length>20||key.length<0) {
                            appUtil.popbox ("序列片段只能输入20个字符，请重新输入!");
                            return false;
                        }
                        return true;
                    }
                }
            }
        }
        //三   基因挖掘搜索
        DataMining.GeneMining = function () {
            var tabNumber = $("#A0253C22EF6447C6ACE32D93D81528BD "+els.tabNumber).length;
            //输入框校验
            var key=$("#annoKeywords").val();
            if(key==""){
                appUtil.popbox ("搜索内容不能为空，请输入后提交!");
                return;
            }
            //功能注释提交按钮
                if (!this.getPerm($("#annoSearchKeyword").attr("permDataId"))) {
                    appUtil.popbox ("您没有使用的权限,请提升用户等级!");
                    return;
                }
                var changeValue = $(els.searchType).val().trim();
                if (changeValue == "SequenceSearch") {//基因序列片段
                    var sequenceCheckBox = $("#sequenceCheckBox").is(":checked"),
                        annoKeywords = $("#annoKeywords").val(),
                        reg = new RegExp("^[TACGtacg]+$", "g");
                    var annoNumber = annoKeywords.split("\n").length;
                    if (tabNumber + annoNumber > _maxTabCount) {
                        appUtil.popbox ("生成的标签过多，不可以提交任务");
                        return;
                    }
                    if (reg.test(annoKeywords)) {
                        $("#annoSearchKeyword").hide();
                        $("#annoSearchKeywordSpan").append("<img src='http://img.biocloud.cn/AppImg/3.1/025.gif'/>");
                        this.postJSONAsync("/report/dge/sequenceSearch", {
                            projectId: _projectId,
                            sequ: annoKeywords,
                            reverse: sequenceCheckBox
                        }, _annoCallbackId(annoKeywords));
                    } else {
                        appUtil.popbox ("输入的序列信息有误,序列中只能包含TACG字母");
                    }
                } else {
                    var keywords = $(els.keywords).val().trim();
                    var type = $(els.searchType).val();
                    if (keywords != "" && _projectId) {
                        if (type === "annoKey") {//功能关键词
                            var keyNumber = keywords.split("\n").length;
                            if (tabNumber + keyNumber > _maxTabCount) {
                                appUtil.popbox ("生成的标签过多，不可以提交任务");
                                return;
                            }
                            $("#annoSearchKeyword").hide();
                            $("#annoSearchKeywordSpan").append("<img src='http://img.biocloud.cn/AppImg/3.1/025.gif'/>");
                            this.postJSONAsync("/report/dge/annoSearchKeyword", {
                                projectId: _projectId,
                                keywords: keywords,
                                type: "annoKey"
                            }, _annoCallback(keywords));
                        } else {
                            //基因编号列表
                            $("#annoSearchKeyword").hide();
                            $("#annoSearchKeywordSpan").append("<img src='http://img.biocloud.cn/AppImg/3.1/025.gif'/>");
                            this.postJSONAsync("/report/dge/annoSearchKeyword", {
                                projectId: _projectId,
                                keywords: keywords,
                                type: "idKey"
                            }, _annoCallbackId(keywords));
                        }
                    } else {
                        appUtil.popbox ("条件不全，不可以提交任务！");
                    }
                }
        };

        //格式化数据
        function _formatData(data) {
            var tr = [];
            for (var i = 0; i < data.length; i++) {
                var td = [];
                td[td.length] = data[i].id;
                var con = data[i].showContent;
                var strings = con.trim().split("\t");
                for (var j = 0; j < strings.length; j++) {
                    var sss = strings[j];
                    td[td.length] = sss;
                }
                tr[tr.length] = td;
            }
            return tr;
        };

        //差异结果查看
        DataMining.degResultLook = function (_this) {
                if (!appUtil.getPerm($(_this).attr("permDataId"))) {
                    appUtil.popbox ("您没有使用的权限,请提升用户等级!");
                    return;
                }
                if(countTab("A5E134F315B144B49F0F6343A304DBA0")){
                    return ;
                }
                var name = $(els.degResultSelect).val();
                dgeselect=$(els.degSelect).val();
                if (name.trim() == "") {
                    appUtil.popbox ("没有具体的差异分组数据!");
                    return;
                }
                $("#degGroupGetDataByType").hide();
                $("#degGroupGetDataByTypeSpan").append("<img src='http://img.biocloud.cn/AppImg/3.1/025.gif'/>");
                var type = "";
                var types = $(els.degResultRegulated);
                for (var i = 0; i < types.length; i++) {
                    if ($(types[i]).is(":checked")==true) {
                        type = $(types[i]).val();
                    }
                }
                if (_rootPath) {

                    DataMining.postJSONAsync("/report/dge/degAnalysis", {
                        path: _rootPath,
                        type: type,
                        dirName: name
                    }, _degCallback(name, type));
                } else {
                    DataMining.warnTipModal("条件不全，不可以提交任务！");
                }
        };

        //差异结果查看回调函数
        function _degCallback(name, type) {
            return function (data) {
                if (data.data.length > 0) {
                    var gUID = _addTabToElement("", "差异分析", 2, false, false, true, 'A5E134F315B144B49F0F6343A304DBA0');
                    $('#data' + gUID).html("");
                    data.data = _formatDegData(data.data);
                    $(els.table).tmpl(data).appendTo('#A5E134F315B144B49F0F6343A304DBA0 #data' + gUID);
                    var map1 = _rootPath + "/DEG_Analysis/" + name + "/DEG_Cluster/" + name + ".png";
                    var map2 = _rootPath + "/DEG_Analysis/" + name + "/" + name + ".cor.png";
                    var imageList=getBase64Images([map1,map2]);
                    var tu = {
                        map1: imageList[0],
                        name1: name + "heatMap.png",
                        name2: name + "sampleMap.png",
                        map2: imageList[1],
                        gUID: gUID,
                        name: name
                    };
                    var jq = $('#A5E134F315B144B49F0F6343A304DBA0 #data' + gUID).parents(".searchResult").find(".mining_mapping_cog").find(".tab-content");
                    $("#headMapAndSoOn").tmpl(tu).appendTo(jq);
                    $("#A5E134F315B144B49F0F6343A304DBA0 .searchResult[relation='" + gUID + "']").find(".pageBtn").attr("fileName", data.fileName);
                    $("#A5E134F315B144B49F0F6343A304DBA0 .searchResult[relation='" + gUID + "']").find("div.js-dataMining-ul-download").attr("fileName", data.fileName);
                    $("#A5E134F315B144B49F0F6343A304DBA0 .searchResult[relation='" + gUID + "']").find(".pageContent").text(data.page);
                    //把结果放到排序筛选中，等待二次搜索
                    $("#A5E134F315B144B49F0F6343A304DBA0 .js-dataMining-sort-analysis").attr("js-deg-type", type);
                    $("#A5E134F315B144B49F0F6343A304DBA0 .searchResult[relation='" + gUID + "']").find(".app_experience_directory #title1").text("差异基因挖掘");
                    $("#A5E134F315B144B49F0F6343A304DBA0 .searchResult[relation='" + gUID + "']").find(".app_experience_directory #title2").text("差异表达基因集查询");
                    if(type==""){
                        type="所有"
                    }else if(type=="up"){
                        type="只上调"
                    }else if(type=="down"){
                        type="只下调"
                    }
                    $("#A5E134F315B144B49F0F6343A304DBA0 .searchResult[relation='" + gUID + "']").find(".app_experience_directory #title3").text("差异分组:"+name+";  类型:"+type);
                    _showContentStyle('A5E134F315B144B49F0F6343A304DBA0');
                } else {
                    DataMining.warnTipModal("没有搜索到结果，请重试！");
                }
                $("#degGroupGetDataByTypeSpan").find("img").hide();
                $("#degGroupGetDataByType").show();
            };
        };
        //获取Base64Image图片,用于一组图片加密
        function getBase64Images(dataList){
            var targetData;
            $.ajax({
                url: "/report/getBase64Image",
                async: false,//改为同步方式
                type: "POST",
                data:{data:JSON.stringify(moreImage(dataList))} ,
                success:function(data) {
                    targetData= data;
                }
            });
            return targetData;
        }
        //差异结果图片base64转码
        function _getImageList(name){
            if(name=="default"){
                //样本图片路径
                var imagePath="http://img.biocloud.cn/cloud_beta_1.0/v_1.3/report/zuoturesult/";
                return [
                    imagePath+"1handle.png",
                    imagePath+"2handle.png",
                    imagePath+"3handle.png",
                    imagePath+"4handle.png",
                    imagePath+"5c.jpg",
                    imagePath+"7handle.png",
                    imagePath+"1.png",
                    imagePath+"1.png"
                ];
            }else {
                //KOG和 KEGGTL是样本图（KEGGTL是一组图片）
                var imageHead=_rootPath+"/DEG_Analysis/"+name;
                var images={
                    "COG":imageHead+"/Cog_Anno/"+name+".Cog.classfy.png",
                    "KOG":"http://img.biocloud.cn/cloud_beta_1.0/v_1.3/report/zuoturesult/2.png",
                    "GO":imageHead+"/go_enrichment/"+name+".GO.png",
                    "KEGGFJ":imageHead+"/Graph/"+name+".KEGG.Phase.png",
                    "GOFJ":imageHead+"/go_enrichment/"+name+".GO.png",//GO富集图
                    "KEGGTL":imageHead+"/pathway/kegg_map/ko00052.png",//kegg通路图
                    "BDL":imageHead+"/DEG_Cluster/hierarchical/"+name+".DEG.cluster.png",//表达量聚类图
                    "YPXG":imageHead+"/"+name+".DEG_cor.png"//样品相关性图
                }
                console.log(images.COG,images.KOG,images.GO,images.KEGGFJ,images.GOFJ,images.KEGGTL,images.BDL,images.YPXG);
                return getBase64Images([images.COG,images.KOG,images.GO,images.KEGGFJ,images.GOFJ,images.KEGGTL,images.BDL,images.YPXG]);
            }
        }
        //遍历json将其转换为jsonobject
        function moreImage(imgPath) {
            var dataList1 = {};
            for (var i = 0; i < imgPath.length; i++) {
                dataList1["imagesrc" + i] = imgPath[i];
            }
            return dataList1
        }
        //差异结果查看数据格式化
        function _formatDegData(data) {
            var tr = [];
            for (var i = 0; i < data.length; i++) {
                var td = [];
                var con = data[i].con;
                var strings = con.trim().split("\t");
                for (var j = 0; j < strings.length; j++) {
                    var sss = strings[j];
                    td[td.length] = sss;
                }
                tr[tr.length] = td;
            }
            return tr;
        };
        //表达基因样品Venn
        DataMining.sampleSelection = function () {
            var samples = Process.samples, js1={};
            js1.items = samples;
            if(!js1.items){
                appUtil.popbox ("没有初始化数据!");
            }
            $("#sequenceSearchPopHtml").tmpl(js1).prependTo("body");
        }
        //表达量挖掘样品信息确认
        DataMining.clickOkSearchFinish = function () {
            geneName = "";
            $("#dajwodaowhdiaskid").find("li").each(function () {
                if ($(this).find("input").is(':checked') == true) {
                    geneName += $(this).find("input").val() + ",";
                }
            });
            if (geneName.substring(0, geneName.length - 1).split(",").length < 2) {
                $("#errorMessage").html("<font color='red'>最少选择两种样品进行比较</font>");
            } else {
                PageStyle.callCode.closePop();
            }
        }
        //表达基因样品限制范围
        /*DataMining.chickGeneSample = function (_this) {
            var dajwodaowhdiaskid = 0;
            $("#dajwodaowhdiaskid li").each(function () {
                if ($(this).find("input").is(':checked') == true) {
                    dajwodaowhdiaskid++;
                }
            });
            if (dajwodaowhdiaskid > 4) {
                $("#dajwodaowhdiaskid li").find("input[checked!='checked']").each(function () {
                    $(this).attr("disabled", "disabled");
                });
            } else {
                $("#dajwodaowhdiaskid li").find("input[checked!='checked']").each(function () {
                    $(this).removeAttr("disabled");
                });
            }
        }*/
        DataMining.chickGeneSample = function () {
            var dajwodaowhdiaskid = 0;
            $("#dajwodaowhdiaskid li").each(function () {
                if ($(this).find("input").is(':checked') == true) {
                    dajwodaowhdiaskid++;
                }
            });
            var check=$("#dajwodaowhdiaskid li").find("input");
            if(dajwodaowhdiaskid > 4){
                for(var i=0;i<check.length;i++){
                    if(($(check[i])).is(":checked")){
                        ($(check[i])).attr("disabled", false);
                    }else{
                        ( $(check[i])).attr("disabled", true);
                    }
                }
            }else{
                $(check).removeAttr("disabled");
            }
        }

        //提交表达量挖掘基因样品和选择信息
        DataMining.sequenceSubmit = function (_this) {
            if (!appUtil.getPerm($(_this).attr("permDataId"))) {
                appUtil.popbox ("您没有使用的权限,请提升用户等级!");
                return;
            }
            //判断是不是有六个标签
            if(countTab("E126A95EB98446739221F07F540FA8AE")){
                return ;
            }
            var geneNum = $("#djaohishdfinloisfsef").val(),
                geneNamepost = geneName.substring(0, geneName.length - 1);
            var items = geneNamepost.split(",");
            if (items.length > 5 || items.length < 2) {
                appUtil.popbox ("选择样品个数为[2,5]");
            } else {
                $("#sequenceSubmit").hide();
                $("#sequenceSubmitSpan").append("<img src='http://img.biocloud.cn/AppImg/3.1/025.gif'/>");
                var gUID = _addTabToElement("", "表达基因", 3, false, true, false, 'E126A95EB98446739221F07F540FA8AE', true);
                $("#E126A95EB98446739221F07F540FA8AE .searchResult[relation='" + gUID + "'] .app_experience_flip").hide();//隐藏分页
                this.postJSONAsync("/report/report/getVennMapData", {
                    path:_rootPath+"/All_gene_expression.xls",
                    geneName: geneNamepost,
                    geneValue: geneNum,
                    degValue:""
                }, DegGroupSelect._vennCallback(geneNum,geneNamepost,gUID, geneName, 'E126A95EB98446739221F07F540FA8AE'));
            }
        }
        //差异分组选择中的差异分组输入框校验
        DataMining.degChangecheck = function (_this) {
            var input = $(_this).val().trim() == null ? $(this).text().trim() : $(_this).val().trim();
            if (input!="") {
                var reg = new RegExp('^[A-Za-z0-9_]+$', 'g');
                if (reg.test(input)) {
                    if (input.length > 4) {
                        appUtil.popbox ("不能超过4个字符的组合");
                        return;
                    }
                } else {
                    appUtil.popbox ("只能输入字母和数字还有下划线");
                    return;
                }
            }else{
                appUtil.popbox ("差异分组样品名称不能为空，请重新输入！");
                return;
            }

        }
        //用来获取注释信息
        DataMining.getAnnotationById = function (id, _this) {
            this.postJSONAsync("/report/dge/report/annotation", {
                id: id,
                projectId: _projectId
            }, _annoById(_this));
        };

        //2基因挖掘
        DataMining.showGeneMining = function (dataGM, gUID) {
            if (dataGM == "" || dataGM == null) {
                $("#A0253C22EF6447C6ACE32D93D81528BD").css("display", "block");
                var gUID = _addTabToElement("", "关键词搜索", 1, false, true, false, 'A0253C22EF6447C6ACE32D93D81528BD');
                $("#A0253C22EF6447C6ACE32D93D81528BD .searchResult[relation='" + gUID + "']").attr("type", "annoKey");
                $("#A0253C22EF6447C6ACE32D93D81528BD .searchResult[relation='" + gUID + "'] .download ").css("display", "none");;
                return gUID;
            } else if (dataGM != null) {
                var url1 = "data:image/png;base64,"+DataQz.data[dataGM.data.length-1].venn;
                var url2 ="data:image/png;base64,"+DataQz.data[dataGM.data.length-1].snv_stat;
                $("#7337140059B34B3BA80698C82488299E .searchResult[relation='" + gUID + "'] .download ").css("display", "block");
                $("#TB_TBJYZCX").removeAttr("disabled");
                $("#TB_TBJYZCX").css("background-image"," linear-gradient(to bottom, #fafdff, #e6e6e6)");
                $("#7337140059B34B3BA80698C82488299E .searchResult[relation='" + gUID + "'] .tableImage1").css("display", "none");
                var Data = {"data": DataQz.data, "gUID": gUID, "fileName": DataQz.fileName};
                $("#7337140059B34B3BA80698C82488299E .searchResult[relation='" + gUID + "']").attr("type", "annoKey");
                $("#tableDataShow").tmpl(Data).appendTo('#7337140059B34B3BA80698C82488299E #data' + gUID);
                ExonUtil.orderListByTd(gUID);
                $("#7337140059B34B3BA80698C82488299E .searchResult[relation='" + gUID + "'] .mining_mapping_cog .mapping_right #pie").find("img").attr("src",url2);
                $("#7337140059B34B3BA80698C82488299E .searchResult[relation='" + gUID + "'] .mining_mapping_cog .mapping_right #venn").find("img").attr("src",url1);
            }
        };

        //根据id获取注释和序列信息
        function _annoById(_this) {
            return function (data) {
                var list = data.annotation.split("\t"), changeValue = $(_this).parents(".searchResult").attr("type"), changeName = $(_this).parents(".searchResult").attr("name");
                var arr = [];
                var databaseColor = [
                    ["COG_class", "#fff"],
                    ["COG_class_annotation", "#F2F5F8"],
                    ["GO_annotation", "#EBEEF0"],
                    ["KEGG_annotation", "#E3E7EB"],
                    ["KOG_class", "#D9DDE0"],
                    ["KOG_class_annotation", "#CBD2D5"],
                    ["Pfam_annotation", "#C6CCCE"],
                    ["Swissprot_annotation", "#BCC4C5"],
                    ["nr_annotation", "#ACBABD"]
                ];
                var annoHtml = "<ul>";
                for (var i = 0; i < list.length; i++) {
                    if (list[i] == "") {
                        continue;
                    }
                    var json = {}, anno, reg = new RegExp(changeName, 'gi');
                    if (/*changeValue == "annoKey" &&*/ changeName != null) {
                        anno = list[i];
                        anno = anno.replace(reg, "<span style='color: #E00FF6'>" + changeName + "</span>");
                    }
                    annoHtml += "<li style='background:" + databaseColor[i][1] + "'><p class=''>" + databaseColor[i][0] + "</p><p>" + (anno == null ? (list[i] == '--' ? (json.con = 'none') : list[i]) : anno) + "</p></li>";
                }
                annoHtml +="</ul>";
                var seq = data.seq;
                if (changeValue == "SequenceSearch" && changeName != null) {
                    var reg = new RegExp(changeName.toUpperCase(), 'g');
                    seq = data.seq.replace(reg, "<span style='color: #E00FF6'>" + changeName.toUpperCase() + "</span>");
                }
                $(_this).parents(".searchResult").find(".right_one").html("");
                $(_this).parents(".searchResult").find(".right_two").html("");
                $(_this).parents(".searchResult").find(".right_two").html(annoHtml);
                $(_this).parents(".searchResult").find(".right_one").html(seq);
            };
        };
        //表格数据操作
        var data_json = {};
        var data_analysis_type = "";
        //序列下载
        DataMining.seqDownload = function () {
            $(els.seqDown).click(function () {
                var fileName = $(this).parents("div").attr("fileName");
                if (!fileName) {
                    appUtil.popbox ("没有可提交的数据,请选择数据再次下载");
                }else{
                    window.open("/report/dge/dataOperationNew?projectId=" + _projectId + "&fileName=" + fileName + "&way=seq");
                }
            });
        };
        DataMining.seqDownload2=function(_this){
            var fileName = $(_this).parents("div").attr("fileName");
            if (!fileName) {
                appUtil.popbox ("没有可提交的数据,请选择数据再次下载");
            }else{
                window.open("/report/dge/dataOperationNew?projectId=" + _projectId + "&fileName=" + fileName + "&way=seq");
            }
        }
        //列表下载
        DataMining.listDownload = function () {
            $(els.listAna).click(function () {
                var fileName = $(this).parents("div").attr("fileName");
                if (!fileName) {
                    appUtil.popbox ("没有可提交的数据,请选择数据再次下载");
                    return;
                } else {
                    window.open("/report/dge/dataOperationNew?fileName=" + fileName + "&projectId=" + _projectId + "&way=list");
                }
            });
        };
        DataMining.listDownload2 = function (_this) {
                var fileName = $(_this).parents("div").attr("fileName");
                if (!fileName) {
                    appUtil.popbox ("没有可提交的数据,请选择数据再次下载");
                    return;
                } else {
                    window.open("/report/dge/dataOperationNew?fileName=" + fileName + "&projectId=" + _projectId + "&way=list");
                }
        };

        DataMining.listAna = function () {
            $(els.listDown).click(function () {
                data_json = _getIdList(this);
                data_analysis_type = 'listAna';
                var html = $(els.analysis).html();
                $(document.body).prepend(html);
            });
        };
        //排序筛选
        DataMining.listAnalysis = function (_this) {
            $(document.body).prepend($("#dataMiningSortAnalysisHTML").html());
            $(".original_pop_up_foot").find("input").val($(_this).attr("js-deg-value"));
        }
        //排序筛选提交
        DataMining.sortAnalysisSubmit = function () {
            var sort = $(".differencesSort li").find("input"),
                manual = $(".differencesManual li"),
                postSort = "", postManual = [],
                change = $(".original_pop_up_foot").find("input").val(),
                type = $(".js-dataMining-sort-analysis").attr("js-deg-type") == null ? "" : $(".js-dataMining-sort-analysis").attr("js-deg-type");

            //排序筛选校验
            if($("#fdr_id").is(":checked")){
                var reg=/^0\.[0-9]{1,5}?$/;/*/^\d+(?=\.{0,1}\d+$|$)/;*/
                var str1=$("#fdr_value").val();
                if( !reg.test(str1)){
                    appUtil.popbox ("您输入的'FDR的值'不正确,请您调整后再提交任务.(注:FDR的范围是(0,0.99999])");
                    return;
                }
            }
            if($("#fc_id").is(":checked")){
                var str1=$("#fc_value").val();
                var reg = new RegExp("^[0-9]+$");
                if( !reg.test(str1)){
                    appUtil.popbox ("您输入的'差异筛选倍数阈值'不正确,请您调整后再提交任务。(注:差异筛选倍数阈值的范围是 >0)");
                    return false;
                }else if (str1 <= 0) {
                    appUtil.popbox ("请输入正确的FC数值必须大于0");
                    return "";
                }
            }

            $(sort).each(function () {
                if ($(this).is(":checked") == true) {
                    postSort = $(this).val();
                }
            });
            if (manual != null) {
                $(manual).each(function () {
                    if ($(this).find("input[name='Experiment_check']").is(":checked") == true) {
                        var checkbox = $(this).find("input[name='Experiment_check']").val();
                        var text = $(this).find("input[name!='Experiment_check']").val();
                        var json = "{'key':'" + checkbox + "','value':'" + text + "'}";
                        postManual.push(json);
                    }
                });
            }
            if(postSort==""&&postManual.join(",")==""){
                PageStyle.callCode.closePop();
                return;
            }
            if (postSort.length > 0 || postManual.length > 0) {
                var json = appUtil.postJSON("/report/dge/sortAnalysis", {
                    projectId: _projectId,
                    sort: postSort,
                    manual: postManual.join(","),
                    change: change,
                    type: type
                });
                PageStyle.callCode.closePop();
                if (json != null) {
                    if (json.error) {
                        var gUID = _addTabToElement("", "排序筛选", 2, false, false, false, 'A5E134F315B144B49F0F6343A304DBA0');
                        _getSortAnalysis(json, gUID, change, type,postSort, 'A5E134F315B144B49F0F6343A304DBA0');
                    }
                }
            } else {
                $(".errorMsg").html("<font color='red'>请选择排序筛选方式</font>");
                //appUtil.popbox ("请选择排序筛选方式");
            }

        }
        //排序筛选回调函数
        function _getSortAnalysis(data, gUID, name, type,sort, htmluuid) {

            if (data.data.length > 0) {
                data.data = _formatDegData(data.data);
                $('#data' + gUID).html("");
                $(els.table).tmpl(data).appendTo('#' + htmluuid + ' #data' + gUID);
                var map1 = _rootPath + "/DEG_Analysis/" + name + "/DEG_Cluster/" + name + ".png";
                var map2 = _rootPath + "/DEG_Analysis/" + name + "/" + name + ".cor.png";
                var tu = {
                    map1: map1,
                    name1: name + "heatMap.png",
                    name2: name + "sampleMap.png",
                    map2: map2,
                    gUID: gUID
                };
                var jq = $('#' + htmluuid + ' #data' + gUID).parents(".searchResult").find(".mining_mapping_cog").find(".tab-content");
                $("#" + htmluuid + " .searchResult[relation='" + gUID + "']").find(".pageBtn").attr("fileName", data.fileName);
                $("#" + htmluuid + " .searchResult[relation='" + gUID + "']").find("div.js-dataMining-ul-download").attr("fileName", data.fileName);
                $("#" + htmluuid + " .searchResult[relation='" + gUID + "']").find(".pageContent").text(data.page);
                $("#headMapAndSoOn").tmpl(tu).appendTo(jq);
                //把结果放到排序筛选中，等待二次搜索
                $("#" + htmluuid + " .js-dataMining-sort-analysis").attr("js-deg-type", type);
                $("#" + htmluuid + " .searchResult[relation='" + gUID + "']").find(".app_experience_directory #title1").text("差异基因挖掘");
                $("#" + htmluuid + " .searchResult[relation='" + gUID + "']").find(".app_experience_directory #title2").text("差异表达基因集查询");
                $("#" + htmluuid + " .searchResult[relation='" + gUID + "']").find(".app_experience_directory #title3").text("排序筛选/ 差异分组:"+name+" 排序: "+sort+" ; ");

                _showContentStyle(htmluuid);
            } else {
                DataMining.warnTipModal("没有搜索到结果，请重试！");
            }
        }

        DataMining.toAnasisly = function () {
            var saveName = $("#data-analysis-name").val();
            var fileId = $("#data-analysis-id").attr("data-pid");
            if (saveName === "" || fileId === "") {
                appUtil.popbox ("文件名和文件路径不可以为空！");
                return;
            }
            window.open("/report/dge/dataOperation?type=" + data_analysis_type + "&json=" + JSON.stringify(data_json) + "&rootPath=" + _rootPath + "&saveName=" + saveName + "&fileId=" + fileId, "_blank");
            PageStyle.callCode.closePop();
        };

        //二次搜索
        DataMining.secondarySearch = function (_this) {
            var type = $(_this).attr("js-data-search-type"), degValue,
                fileName = $(_this).parents(".searchResult").find("div.js-dataMining-ul-download").attr("filename"),
                keyword = $(_this).parents(".searchResult").find(".js-data-table-search").val().trim(),
                relation=$(_this).parents(".searchResult").attr("Relation");
            //二次搜索输入框校验
            if(keyword==""){
                appUtil.popbox ("请输入需要查询的注释信息");
                return;
            }else{
                var reg=new RegExp("^[A-Za-z0-9]+$");
                if(!reg.test(keyword)){
                    appUtil.popbox ("搜索只能有字母和数字，请重新输入!");
                    return;
                } else {
                    if(keyword.length>20||keyword.length<0) {
                        appUtil.popbox ("搜索只能输入20个字符，请重新输入!");
                        return;
                    }
                }
            }
            if (fileName.trim() == "") {
                appUtil.popbox ("没有基因数据!");
                return;
            }
            $(_this).siblings(".sousuo"+relation).append("<img src='http://img.biocloud.cn/AppImg/3.1/025.gif'/>");
            $(_this).hide();
            this.postJSONAsync("/report/dge/AllSecondarySearch", {
                projectId: _projectId,
                keyword: keyword,
                domainfileName: fileName
            }, _secondSearch(_this));
        };
        //二次搜索回调函数
        function _secondSearch(_this) {
            return function (data) {
                if (data) {
                    var tab = $(_this).parents(".searchResult").find(".mining_tbale"),
                        gUID = $(_this).parents(".searchResult").attr("Relation"),
                        keywords = $(_this).parents(".searchResult").find(".js-data-table-search").val();
                    if (data.data.length > 1) {
                        data.data = _formatDegData(data.data);
                        tab.html("");
                        $(_this).parents(".searchResult[relation='" + gUID + "']").attr("name", keywords);
                        $(_this).parents(".searchResult").find("div.js-dataMining-ul-download").attr("filename",data.fileName);
                        $(_this).parents(".searchResult[relation='" + gUID + "']").find(".pageContent").text(data.page);
                        $(els.table).tmpl(data).appendTo(tab);
                    } else {
                        appUtil.popbox ("没有搜索到结果，请重试！");
                    }
                    $(_this).siblings(".sousuo"+gUID).find("img").hide();
                    $(_this).show();
                }
            };
        };

        //用来获取表格中用户选中的id列表
        function _getIdList(_this) {
            var table = $(_this).parents(".searchResult").find("table");
            var fileName = table.attr("filename");
            var json = {isAll: false, list: [], fileName: fileName};
            if (table.find("input[name=selectAll]").attr("checked")) {
                json.isAll = true;
            } else {
                var ids = table.find("input[name=selectOne]");
                for (var int = 0, l = ids.length; int < l; int++) {
                    if ($(ids[int]).attr("checked")) {
                        json.list[json.list.length] = $(ids[int]).val();
                    }
                }
            }
            return json;
        };

        //全选按钮的操作
        DataMining.clickSelectAll = function (_this) {
            if ($(_this).attr("checked")) {
                $(_this).parents("table").find("input[type=checkbox]").attr("checked", true);
            } else {
                $(_this).parents("table").find("input[type=checkbox]").attr("checked", false);
            }
        };

        //用来显示大图片
        DataMining.showImageBig = function (_this) {
            var str=$(_this).find("img").attr("src");
            if(str=="http://img.biocloud.cn/cloud_beta_1.0/v_1.3/report/report_loading.gif"){

            }else{
                var gUID = $(_this).parents(".searchResult").attr("markAttr");
                var path = $(_this).find("img").attr("src");
                var selectattr = $(_this).find("img").attr("selectattr");
                if (selectattr == "heatMap") {
                    var data = {path: path, gUID: gUID};
                    $("#showBigPngImage").tmpl(data).prependTo("body");
                    return;
                }
                var path11 = path.replace(".png", ".svg");
                path11 = (path11.split("=")[1] ? path11.split("=")[1] : path11.split("=")[0]);
                var flag = this.postJSON("/report/dge/isCompletePic", {path: path11});
                if (flag) {
                    path = path.replace(".svg", ".png");
                    var data = {path: path, gUID: gUID};
                    $("#showBigPngImage").tmpl(data).prependTo("body");
                } else {
                    var data = {path: path, gUID: gUID};
                    $("#showBigPngImage").tmpl(data).prependTo("body");
                }
            }
        };

        //移除大图
        DataMining.hideImageBig = function () {
            $(".tab_big").remove();
            $(".modal").remove();
            $(".modal-backdrop").hide();
        };

        //svg 图片点击触发代码
        DataMining.vennClick = function (path, keyword, _type) {
            var gUID = $(".tab_big").attr("markAttr");
            this.postJSONAsync("/report/dge/svgHeatClick", {
                rootPath: _rootPath,
                path: path,
                keyword: keyword,
                type: _type
            }, _svgImageClickCallback(gUID));
        };

        //svg图片点击热区，获取数据的回调函数
        function _svgImageClickCallback(gUID) {
            return function (data) {
                if (data.data.length > 0) {

                    data.data = _formatData(data.data)
                    var asdf = $('.searchResult[markAttr=' + gUID + ']').find("#data" + gUID);
                    asdf.html("");
                    $("#tableDataShow").tmpl(data).appendTo(asdf);
                    DataMining.hideImageBig();
                } else {
                    DataMining.warnTipModal("没有搜索到结果，请重试！");
                }
            };
        };

        //cog kog 图片热点点击
        DataMining.cogKogClick = function (path, num, keyword) {
            var gUID = $(".tab_big").attr("markAttr");
            this.postJSONAsync("/report/dge/cogOrKogSvgHeatClick", {
                path: path,
                keyword: keyword,
                num: num
            }, _cogKogSvgImageClickCallback(gUID));
        };

        function _cogKogSvgImageClickCallback(gUID) {
            return function (data) {
                if (data.length > 0) {

                    var s = new Object();
                    s.data = _format(data);
                    var asdf = $('.searchResult[markAttr=' + gUID + ']').find("#data" + gUID);
                    asdf.html("");
                    $("#tableDataShow").tmpl(s).appendTo(asdf);
                    DataMining.hideImageBig();
                } else {
                    appUtil.popbox ("没有搜索到结果，请重试！");
                }
            };
        };

        function _format(data) {
            var tr = [];
            for (var i = 0; i < data.length; i++) {
                var td = [];
                var con = data[i].con;
                var strings = con.trim().split("\t");
                for (var j = 0; j < strings.length; j++) {
                    var sss = strings[j];
                    td[td.length] = sss;
                }
                tr[tr.length] = td;
            }
            return tr;
        };

        var DegGroupSelect = (function () {
            var DegGroupSelect = new Object();
            _extends(DegGroupSelect, Util);
            var elements = {
                "popHtml": "#degCombinationSelectPopHtml", //弹出框html模板
                "popClick": "#degCombinationSelectPop", //差异结果分析弹出框按钮id
                "input": "#dee39c input[type=checkbox]"
            };

            //初始化差异组合筛选弹出框
            DegGroupSelect.init = function () {
                if (DataMining.degGroupSelectArr.length > 0) {
                    var js = {};
                    var items = DataMining.degGroupSelectArr;
                    js.items = items;
                    $(elements.popHtml).tmpl(js).prependTo("body");
                    //用来保证用户只能选择5个表达式
                    (function () {
                        $("#dee39c").delegate("input[type=checkbox]", "click", function () {
                            var selectVennys = $("#dee39c").find("input[type=checkbox]");
                            var checkedSelects = [];
                            var nocheckedSelects = [];
                            for (var int = 0; int < selectVennys.length; int++) {
                                var qweqwe = selectVennys[int];
                                (qweqwe.checked) ? (checkedSelects[checkedSelects.length] = qweqwe) : (nocheckedSelects[nocheckedSelects.length] = qweqwe);
                            }
                            if (checkedSelects.length == 5) {
                                for (var int = 0, l = nocheckedSelects.length; int < l; int++) {
                                    var nocheckedSelect = nocheckedSelects[int];
                                    $(nocheckedSelect).attr("disabled", true);
                                }
                            } else {
                                for (var int = 0, l = nocheckedSelects.length; int < l; int++) {
                                    var nocheckedSelect = nocheckedSelects[int];
                                    $(nocheckedSelect).attr("disabled", false);
                                }
                            }
                        });
                    })();
                } else {
                    appUtil.popbox ("该分析只对两种或以上差异组合生效");
                }
            };

            //差异组合选择
            //第四个按钮下面的 差异组合选择按钮
            DegGroupSelect.openPop = function () {
                $("#degCombinationSelectPop").click(function () {
                    if (!_isBlank) {
                        DegGroupSelect.init();
                    } else {
                        DegGroupSelect.warnTipModal("您还没有选择结题报告，无法进行其他操作");
                    }
                });
            };

            //差异组合最终选择
            DegGroupSelect.clickOkDegCombinationFinish = function () {
                var jQueryItem = $("#dee39c");
                var items = jQueryItem.find("input[type=checkbox]");
                var checkedSelects = [];
                for (var i = 0; i < items.length; i++) {
                    if ($(items[i]).is(':checked')) {
                        checkedSelects[checkedSelects.length] = items[i];
                    }
                }
                var combination = $("#dee39c").find(" li"), comName = "",num=0;
                $(combination).each(function () {
                    if ($(this).find("input[type='checkbox']").is(":checked")==true) {
                        if(($(this).find(".checkbox").find("label").find(".differencesCom").val())==""){
                            appUtil.popbox("差异组合样品名称不能为空");
                            error=true;
                            num=1;
                            return;
                        }else if(($(this).find(".checkbox").find("label").find(".differencesCom").val()).length>4){
                            appUtil.popbox("差异组合样品名称不能超过4个字符");
                            error=true;
                            num=1;
                            return;
                        }
                        comName += $(this).find(".checkbox").find("label").find(".differencesCom").val() + ",";
                    }
                });
                this.degCombinationFinishCheckedSelects = checkedSelects;
                this.comName = comName;
                if(num==0){
                    PageStyle.callCode.closePop();
                }
            };

            //差异结果分析提交按钮触发代码
            DegGroupSelect.submit = function () {
                $("#eddc65964cdc8e49fcb59433dcf1").click(function () {
                    if (!appUtil.getPerm($(this).attr("permDataId"))) {
                        appUtil.popbox ("您没有使用的权限,请提升用户等级!");
                        return;
                    }
                    if(countTab("A5E134F315B144B49F0F6343A304DBA0")){
                        return ;
                    }
                    var cfg = {};
                    var checkedSelects = DegGroupSelect.degCombinationFinishCheckedSelects;
                    var gUID = null;
                    if(error){
                        appUtil.popbox ("参数不对，不能提交!");
                        return;
                    }else{
                        if (checkedSelects && checkedSelects.length > 1) {
                            $("#eddc65964cdc8e49fcb59433dcf1").hide();
                            $("#eddc65964cdc8e49fcb59433dcf1Span").append("<img src='http://img.biocloud.cn/AppImg/3.1/025.gif' />");
                            cfg.venn = checkedSelects;
                            var venns = "";
                            for (var j = 0; j < checkedSelects.length; j++) {
                                if (j == 0) {
                                    venns += $(checkedSelects[j]).val()+"_regulated";
                                } else {
                                    venns += "," + $(checkedSelects[j]).val()+"_regulated";
                                }
                            }
                            var type = "all";
                            var types = $("input[name=degSelect]");
                            for (var i = 0; i < types.length; i++) {
                                if ($(types[i]).is(":checked")==true) {
                                    type = $(types[i]).val();
                                }
                            }
                            cfg.type = type;
                            gUID = _addTabToElement("", "差异分组", 3, false, false, false, 'A5E134F315B144B49F0F6343A304DBA0');
                            $("#A5E134F315B144B49F0F6343A304DBA0 .searchResult[relation='" + gUID + "'] .app_experience_flip").hide();//隐藏分页
                            appUtil.postJSONAsync("/report/report/getVennMapData", {
                                path: _rootPath+"/All_gene_expression.xls",
                                geneName: venns,
                                geneValue:"",
                                degValue: type
                            }, DegGroupSelect._vennCallback(venns,type,gUID, DegGroupSelect.comName.substring(0, DegGroupSelect.comName.length - 1), 'A5E134F315B144B49F0F6343A304DBA0'));
                        } else {
                            appUtil.popbox ("该分析只对两种或以上差异组合生效");
                        }
                    }

                });
            };
            //这里是venn图的回调函数
            DegGroupSelect._vennCallback = function (geneNum,type,gUID, comName, htmluuid) {

                return function(data){
                    if (data === null){
                        appUtil.popbox ("Sorry！您的维恩图做图制作过程出现问题，我们会努力解决中，请谅解。");
                        return;
                    }else{
                        appUtil.popbox ("venn图任务提交成功，请耐心等候几分钟。");
                        setTimeout(DegGroupSelect.isVennSuccess(data.resultPath,gUID,comName,htmluuid,geneNum,type), 5000);
                    }
                };
            };

            //判断Venn图是否完成
            DegGroupSelect.isVennSuccess = function (resultPath,gUID,comName,htmluuid,geneNum,type){
                appUtil.postJSONAsync(_projectName+"/report/report/isVennComplete",{path:resultPath}, _isVennFinish(resultPath,gUID,comName,htmluuid,geneNum,type));
            }
            //判断Venn图是否完成回调函数
            function _isVennFinish(resultPath,gUID,comName,htmluuid,geneNum,type){
                return function(data){
                    if(data.toString()=="false" ){
                        setTimeout(function(){
                            DegGroupSelect.isVennSuccess(resultPath,gUID,comName,htmluuid,geneNum,type)
                        }, 5000);
                    }else {
                        DegGroupSelect.uuid = data.uuid;
                        var size = data.size;
                        switch (size) {
                            case 3:
                                _showVenn(gUID, 0)
                                $(".Map1_1" + gUID).find("span").text(data.map_1);
                                $(".Map1_2" + gUID).find("span").text(data.map_2);
                                $(".Map1_all" + gUID).find("span").text(data.map_all);
                                break;
                            case 7:
                                _showVenn(gUID, 1)
                                $(".Map2_1" + gUID).find("span").text(data.map_1);
                                $(".Map2_2" + gUID).find("span").text(data.map_2);
                                $(".Map2_3" + gUID).find("span").text(data.map_3);
                                $(".Map2_all" + gUID).find("span").text(data.map_all);
                                $(".Map2_12" + gUID).find("span").text(data.map_12);
                                $(".Map2_13" + gUID).find("span").text(data.map_13);
                                $(".Map2_23" + gUID).find("span").text(data.map_23);
                                break;
                            case 15:
                                _showVenn(gUID, 2)
                                $(".Map3_1" + gUID).find("span").text(data.map_1);
                                $(".Map3_2" + gUID).find("span").text(data.map_2);
                                $(".Map3_3" + gUID).find("span").text(data.map_3);
                                $(".Map3_4" + gUID).find("span").text(data.map_4);
                                $(".Map3_all" + gUID).find("span").text(data.map_all);
                                $(".Map3_12" + gUID).find("span").text(data.map_12);
                                $(".Map3_13" + gUID).find("span").text(data.map_13);
                                $(".Map3_14" + gUID).find("span").text(data.map_14);
                                $(".Map3_23" + gUID).find("span").text(data.map_23);
                                $(".Map3_24" + gUID).find("span").text(data.map_24);
                                $(".Map3_34" + gUID).find("span").text(data.map_34);
                                $(".Map3_123" + gUID).find("span").text(data.map_123);
                                $(".Map3_124" + gUID).find("span").text(data.map_124);
                                $(".Map3_134" + gUID).find("span").text(data.map_134);
                                $(".Map3_234" + gUID).find("span").text(data.map_234);
                                break;
                            case 31:
                                _showVenn(gUID, 3)
                                $(".Map4_1" + gUID).find("span").text(data.map_1);
                                $(".Map4_2" + gUID).find("span").text(data.map_2);
                                $(".Map4_3" + gUID).find("span").text(data.map_3);
                                $(".Map4_4" + gUID).find("span").text(data.map_4);
                                $(".Map4_5" + gUID).find("span").text(data.map_5);
                                $(".Map4_all" + gUID).find("span").text(data.map_all);
                                $(".Map4_12" + gUID).find("span").text(data.map_12);
                                $(".Map4_13" + gUID).find("span").text(data.map_13);
                                $(".Map4_14" + gUID).find("span").text(data.map_14);
                                $(".Map4_15" + gUID).find("span").text(data.map_15);
                                $(".Map4_23" + gUID).find("span").text(data.map_23);
                                $(".Map4_24" + gUID).find("span").text(data.map_24);
                                $(".Map4_25" + gUID).find("span").text(data.map_25);
                                $(".Map4_34" + gUID).find("span").text(data.map_34);
                                $(".Map4_35" + gUID).find("span").text(data.map_35);
                                $(".Map4_45" + gUID).find("span").text(data.map_45);
                                $(".Map4_123" + gUID).find("span").text(data.map_123);
                                $(".Map4_124" + gUID).find("span").text(data.map_124);
                                $(".Map4_125" + gUID).find("span").text(data.map_125);
                                $(".Map4_134" + gUID).find("span").text(data.map_134);
                                $(".Map4_135" + gUID).find("span").text(data.map_135);
                                $(".Map4_145" + gUID).find("span").text(data.map_145);
                                $(".Map4_234" + gUID).find("span").text(data.map_234);
                                $(".Map4_235" + gUID).find("span").text(data.map_235);
                                $(".Map4_245" + gUID).find("span").text(data.map_245);
                                $(".Map4_345" + gUID).find("span").text(data.map_345);
                                $(".Map4_1234" + gUID).find("span").text(data.map_1234);
                                $(".Map4_1235" + gUID).find("span").text(data.map_1235);
                                $(".Map4_1245" + gUID).find("span").text(data.map_1245);
                                $(".Map4_1345" + gUID).find("span").text(data.map_1345);
                                $(".Map4_2345" + gUID).find("span").text(data.map_2345);
                                break;
                            default:
                                break;
                        }
                        _getImageByJson(size, data, gUID, comName, htmluuid);
                        _showContentStyle(htmluuid);
                        if (htmluuid == "E126A95EB98446739221F07F540FA8AE") {
                            $("#E126A95EB98446739221F07F540FA8AE .searchResult[relation='" + gUID + "']").find("#loadImg1").hide();
                            $("#E126A95EB98446739221F07F540FA8AE .searchResult[relation='" + gUID + "']").find("#loadImg").hide();
                            $("#E126A95EB98446739221F07F540FA8AE .searchResult[relation='" + gUID + "']").find(".app_experience_directory #title1").text("表达量挖掘");
                            $("#E126A95EB98446739221F07F540FA8AE .searchResult[relation='" + gUID + "']").find(".app_experience_directory #title2").text("样品表达基因集维恩图");
                            $("#E126A95EB98446739221F07F540FA8AE .searchResult[relation='" + gUID + "']").find(".app_experience_directory #title3").text("参数: "+geneNum+" ; "+type);
                            $("#sequenceSubmitSpan").find("img").hide();
                            $("#sequenceSubmit").show();
                        }
                        if (htmluuid == "A5E134F315B144B49F0F6343A304DBA0") {
                            $("#A5E134F315B144B49F0F6343A304DBA0 .searchResult[relation='" + gUID + "']").find("#loadImg1").hide();
                            $("#A5E134F315B144B49F0F6343A304DBA0 .searchResult[relation='" + gUID + "']").find("#loadImg").hide();
                            //差异表达基因集维恩图
                            $("#A5E134F315B144B49F0F6343A304DBA0 .searchResult[relation='" + gUID + "']").find(".app_experience_directory #title1").text("差异基因挖掘");
                            $("#A5E134F315B144B49F0F6343A304DBA0 .searchResult[relation='" + gUID + "']").find(".app_experience_directory #title2").text("差异表达基因集维恩图");
                            if(type=="all"){
                                type="所有"
                            }else if(type=="up"){
                                type="只上调"
                            }else if(type=="down"){
                                type="只下调"
                            }
                            $("#A5E134F315B144B49F0F6343A304DBA0 .searchResult[relation='" + gUID + "']").find(".app_experience_directory #title3").text("差异分组: "+geneNum+";  类型:"+type);
                            $("#eddc65964cdc8e49fcb59433dcf1Span").find("img").hide();
                            $("#eddc65964cdc8e49fcb59433dcf1").show();
                        }
                    }
                };
            };
            //添加venn图显示数字
            function _getImageByJson(size, json, gUID, comName, htmluuid) {
                var imagePath, dataText = [], numFont, usemap;
                switch (size) {
                    case 3:
                        usemap = "1";
                        imagePath ="/static/js/inside/report/vennimg/2W.png";
                        numFont = '20pt';
                        dataText.push({num: json.map_2, x: 65, y: 122});
                        dataText.push({num: json.map_all, x: 170, y: 122});
                        dataText.push({num: json.map_1, x: 290, y: 122});
                        dataText.push({num: comName.split(",")[0], x: 328, y: 30, rotate: Math.PI / 4});
                        dataText.push({num: comName.split(",")[1], x: -190, y: 215, rotate: -Math.PI / 2});
                        break;
                    case 7:
                        usemap = "2";
                        imagePath ="/static/js/inside/report/vennimg/3W.png";
                        numFont = '15pt';
                        dataText.push({num: json.map_1, x: 242, y: 69});
                        dataText.push({num: json.map_2, x: 159, y: 202});
                        dataText.push({num: json.map_3, x: 82, y: 69});
                        dataText.push({num: json.map_12, x: 218, y: 144});
                        dataText.push({num: json.map_13, x: 159, y: 63});
                        dataText.push({num: json.map_23, x: 117, y: 144});
                        dataText.push({num: json.map_all, x: 167, y: 122});
                        dataText.push({num: comName.split(",")[1], x: 262, y: 217});
                        dataText.push({num: comName.split(",")[0], x: 290, y: 20, rotate: Math.PI / 4});
                        dataText.push({num: comName.split(",")[2], x: -150, y: 170, rotate: -Math.PI / 2});
                        break;
                    case 15:
                        usemap = "3";
                        imagePath ="/static/js/inside/report/vennimg/4W.png";
                        numFont = '10pt';
                        dataText.push({num: json.map_1, x: 285, y: 108});
                        dataText.push({num: json.map_2, x: 212, y: 29});
                        dataText.push({num: json.map_3, x: 126, y: 33});
                        dataText.push({num: json.map_4, x: 68, y: 106});
                        dataText.push({num: json.map_12, x: 245, y: 64});
                        dataText.push({num: json.map_13, x: 248, y: 144});
                        dataText.push({num: json.map_14, x: 170, y: 210});
                        dataText.push({num: json.map_23, x: 168, y: 55});
                        dataText.push({num: json.map_24, x: 101, y: 144});
                        dataText.push({num: json.map_34, x: 100, y: 70});
                        dataText.push({num: json.map_123, x: 212, y: 96});
                        dataText.push({num: json.map_124, x: 135, y: 178});
                        dataText.push({num: json.map_134, x: 205, y: 182});
                        dataText.push({num: json.map_234, x: 123, y: 103});
                        dataText.push({num: json.map_all, x: 166, y: 140});
                        dataText.push({num: comName.split(",")[0], x: 321, y: 81});
                        dataText.push({num: comName.split(",")[1], x: 274, y: 24});
                        dataText.push({num: comName.split(",")[2], x: 71, y: 28});
                        dataText.push({num: comName.split(",")[3], x: 24, y: 81});
                        break;
                    case 31:
                        usemap = "4";
                        imagePath ="/static/js/inside/report/vennimg/5W.png";
                        numFont = '8pt';
                        dataText.push({num: json.map_1, x: 187, y: 27});
                        dataText.push({num: json.map_2, x: 274, y: 114});
                        dataText.push({num: json.map_3, x: 211, y: 210});
                        dataText.push({num: json.map_4, x: 115, y: 192});
                        dataText.push({num: json.map_5, x: 92, y: 81});
                        dataText.push({num: json.map_12, x: 231, y: 72});
                        dataText.push({num: json.map_13, x: 213, y: 180});
                        dataText.push({num: json.map_14, x: 189, y: 55});
                        dataText.push({num: json.map_15, x: 149, y: 52});
                        dataText.push({num: json.map_23, x: 245, y: 156});
                        dataText.push({num: json.map_24, x: 131, y: 167});
                        dataText.push({num: json.map_25, x: 248, y: 107});
                        dataText.push({num: json.map_34, x: 163, y: 193});
                        dataText.push({num: json.map_35, x: 121, y: 96});
                        dataText.push({num: json.map_45, x: 107, y: 129});
                        dataText.push({num: json.map_123, x: 219, y: 147});
                        dataText.push({num: json.map_124, x: 203, y: 67});
                        dataText.push({num: json.map_125, x: 218, y: 94});
                        dataText.push({num: json.map_134, x: 192, y: 174});
                        dataText.push({num: json.map_135, x: 140, y: 86});
                        dataText.push({num: json.map_145, x: 167, y: 76});
                        dataText.push({num: json.map_234, x: 162, y: 165});
                        dataText.push({num: json.map_235, x: 240, y: 126});
                        dataText.push({num: json.map_245, x: 129, y: 145});
                        dataText.push({num: json.map_345, x: 128, y: 117});
                        dataText.push({num: json.map_1234, x: 193, y: 152});
                        dataText.push({num: json.map_1235, x: 215, y: 118});
                        dataText.push({num: json.map_1245, x: 187, y: 88});
                        dataText.push({num: json.map_1345, x: 150, y: 104});
                        dataText.push({num: json.map_2345, x: 152, y: 142});
                        dataText.push({num: json.map_all, x: 180, y: 121});
                        dataText.push({num: comName.split(",")[0], x: 247, y: 15});
                        dataText.push({num: comName.split(",")[1], x: 308, y: 113});
                        dataText.push({num: comName.split(",")[2], x: 274, y: 218});
                        dataText.push({num: comName.split(",")[3], x: 71, y: 218});
                        dataText.push({num: comName.split(",")[4], x: 68, y: 57});
                        break;
                }

                var img = new Image();
                img.onload = function () {
                    var canvas = convertImageToCanvas(img, dataText, numFont);
                    img.src = canvas.toDataURL();
                    $("#" + htmluuid + " .searchResult[Relation='" + gUID + "'] div.imageAddCanvasData").html("<img style='display: inline-block;max-width: none;max-height:none' src='" + img.src + "' usemap='#Map" + usemap + gUID + "''/>");
                    //$("#" + htmluuid + " .searchResult[Relation='" + gUID + "'] div.imageAddCanvasData").html("<img src='" + img.src + "' usemap='#Map" + usemap + gUID + "''/>");
                    img.onload = null;
                };
                img.src = imagePath;
            }

            //选择区域，获取venn图的列表
            DegGroupSelect.showVennyTable = function (_this) {
                if($(_this).next().find("span").text()==0){
                    appUtil.popbox("没有结果!");
                    return ;
                }
                var mark = $(_this).attr("mark");
                this.postJSONAsync("/report/dge/showVennyTable", {
                    uuid: this.uuid,
                    projectId: _projectId,
                    key: mark
                }, _showVennyCallback(_this));
            };

            //Venn复选区域
            DegGroupSelect.changeVennCheck = function (_this) {
                var key = "";
                $(_this).parents(".problem_map").find("div input").each(function () {
                    if ($(this).attr("checked") == "checked") {
                        key += $(this).attr("mark") + ",";
                    }
                });
                if (key != null && key.length > 0) {
                    this.postJSONAsync("/report/dge/showVennyTable", {
                        uuid: this.uuid,
                        projectId: _projectId,
                        key: key.substring(0, key.length - 1)
                    }, _showVennyCallback(_this));
                }
            }

            //显示venn图的回调函数
            function _showVennyCallback(_this) {
                return function (data) {
                    if(JSON.stringify(data)=="{}"){
                        appUtil.popbox ("没有搜索到结果，请重试！");
                    }else{
                        data.data = _formatData(data.data)
                        var gUID = $(_this).parents(".searchResult").attr("markAttr");
                        $('#data' + gUID).html("");
                        $(".searchResult[relation='" + gUID + "']").find(".pageBtn").attr("fileName", data.fileName);
                        $(".searchResult[relation='" + gUID + "']").find("div.js-dataMining-ul-download").attr("fileName", data.fileName);
                        $(".searchResult[relation='" + gUID + "']").find(".pageContent").text(data.page);
                        $(".searchResult[relation='" + gUID + "']").find(".pagenow_text").text("1");
                        $(".searchResult[relation='" + gUID + "'] .app_experience_flip").show();//隐藏分页
                        $("#tableDataShow").tmpl(data).appendTo('#data' + gUID);
                    }
                }
            };

            //展示venn图
            function _showVenn(gUID, number) {
                var lis = $("#vennMap" + gUID).find("a");
                for (var int = 0; int < lis.length; int++) {
                    if (int == number) {
                        var img = $(lis[int]).find("img").attr("src");
                        $(".searchResult[relation='" + gUID + "']").find("#downloadVenninput").attr("onclick", "DgeReport.DataMining.DegGroupSelect.downloadVenn('" + int + "','" + gUID + "')");
                        $(lis[int]).show();
                    } else {
                        $(lis[int]).hide();
                    }
                }
            };
            //给图片加载数据
            function convertImageToCanvas(image, json, numFont) {
                var canvas = document.createElement("canvas");
                canvas.width = image.width;
                canvas.height = image.height;
                canvas.getContext("2d").drawImage(image, 0, 0);
                for (var item in json) {
                    canvas.getContext("2d").font = numFont + " Calibri";  //文字字体
                    var it = json[item];
                    canvas.getContext("2d").save();
                    if (it.rotate) {
                        canvas.getContext("2d").translate(it.x, it.y);
                        canvas.getContext("2d").rotate(it.rotate);
                        if(it.num) {
                            if (it.num.length > 5) {
                                canvas.getContext("2d").font = "6pt Calibri";  //文字字体
                            }
                        }
                        canvas.getContext("2d").fillText(it.num, 0, 0);
                    } else {
                        canvas.getContext("2d").fillText(it.num, it.x, it.y);
                    }
                }
                return canvas;
            }
            /*下载图片*/
            function loadImage(path) {
                var image = new Image();
                image.src = path;
                image.src = image.src.replace("image/png", 'image/octet-stream');
                var saveFile = function (data, filename) {
                    var save_link = document.createElementNS('http://www.w3.org/1999/xhtml', 'a');
                    save_link.download = filename;
                    save_link.href = data;
                    if (isIE()) {
                        var canvas = document.createElement("canvas");
                        canvas.width = image.width;
                        canvas.height = image.height;
                        canvas.getContext("2d").drawImage(image, 0, 0);
                        canvas.toBlob(function (blob) {
                            saveAs(blob, filename);
                        });
                    } else {
                        var event;
                        if (document.createEvent) {
                            event = document.createEvent('MouseEvents');
                            event.initMouseEvent('click', false, false, window, 0, 0, 0, 0, 0, false, false, false, false, 0, null);
                            save_link.dispatchEvent(event);
                        }
                    }
                };
                var filename = 'appUtil_' + (new Date()).getTime() + '.' + "png";
                saveFile(image.src, filename);
            }
            /*下载相关性图*/
            DegGroupSelect.downloadDependency = function (_this, gUID) {
                var path=$(_this).parent("div").next("div.img_data").find("img").attr("src");
                if(path=="http://img.biocloud.cn/cloud_beta_1.0/v_1.3/report/report_loading.gif"){
                    appUtil.popbox ("绘图任务没有完成，不能进行下载!");
                    return;
                }else{
                    loadImage(path);
                }

            }
            //下载venn图片
            DegGroupSelect.downloadVenn = function (num, gUID) {
                var image = new Image();
                image.src = $(".searchResult[Relation='" + gUID + "'] div.imageAddCanvasData").find("img")[0].src;
                image.src = image.src.replace("image/png", 'image/octet-stream');
                var saveFile = function (data, filename) {
                    var save_link = document.createElementNS('http://www.w3.org/1999/xhtml', 'a');
                    save_link.download = filename;
                    save_link.href = data;

                    if (isIE()) {
                        var canvas = document.createElement("canvas");
                        canvas.width = image.width;
                        canvas.height = image.height;
                        canvas.getContext("2d").drawImage(image, 0, 0);
                        canvas.toBlob(function (blob) {
                            saveAs(blob, filename);
                        });
                    } else {
                        var event;
                        if (document.createEvent) {
                            event = document.createEvent('MouseEvents');
                            event.initMouseEvent('click', false, false, window, 0, 0, 0, 0, 0, false, false, false, false, 0, null);
                            save_link.dispatchEvent(event);
                        }
                    }
                };
                var filename = 'appUtil_' + (new Date()).getTime() + '.' + "png";
                saveFile(image.src, filename);
            }
            return DegGroupSelect;
        })();
        //判断是不是该死的IE
        var isIE = function () {
            var Sys = {};
            var ua = navigator.userAgent.toLowerCase();
            var s;
            (s = ua.match(/rv:([\d.]+)\) like gecko/)) ? Sys.ie = s[1] :
                (s = ua.match(/msie ([\d.]+)/)) ? Sys.ie = s[1] : 0;
            if (Sys.ie) {
                return true;
            } else {
                return false;
            }
        }
        //数据绘图
        var DrawMapping = (function () {
            var DrawMapping = new Object();
            _extends(DrawMapping, Util);

            //画图任务提交
            DrawMapping.drawingMaps = function (type, _this) {
                var table = $(_this).parents(".searchResult").find("div.js-dataMining-ul-download");
                var fileName = table.attr("filename");
                if(!fileName){
                    appUtil.popbox ("请选择基因数据!");
                    return;
                }
                var idList = {isAll: true, list: []}, filePath;
                idList.id = this.getUID;
                idList.type = type;
                idList.fileName = fileName;
                idList.rootPath = _rootPath;
                idList.change =$(_this).parents(".searchResult").find(".js-dataMining-sort-analysis").attr("js-deg-value") == null ? "" : $(_this).parents(".searchResult").find(".js-dataMining-sort-analysis").attr("js-deg-value");
                if (_drawMapPath) {
                    idList.path = _drawMapPath + "/" + this.getUID();
                    idList.name = this.getUID();
                } else {
                    throw "没有初始化画图图片缓存路径";
                }
                if (idList.type == 'kegg') {
                    filePath = idList.path + "/pathway/kegg_enrichment/" + idList.name;
                } else if (idList.type == 'go') {
                    filePath = idList.path + "/go_enrichment/" + idList.name ;
                } else if (idList.type == 'deggo') {
                    filePath = idList.path + "/go_enrichment/" + idList.name + ".GO";
                } else if (idList.type == 'tonglu') {
                    idList.dgeselect=dgeselect;
                    filePath = idList.path + "/pathway/kegg_map/";
                } else if (idList.type == 'tgb') {
                    idList.path += "/TMP_OD/";
                    idList.dgeselect=dgeselect;
                    filePath = idList.path + "DEG_Cluster/" +"heatmap";
                } else {
                    filePath = idList.path + "/" + idList.name;
                }
                if (!idList.isAll && idList.list.length == 0) {
                    appUtil.popbox ("你还没有选择数据，无法提交画图任务。")
                    return;
                }
                appUtil.popbox ("画图任务提交成功，请耐心等候几分钟。");
                var str=$(_this).parents(".searchResult").find("img[selectAttr=" + (type == "tgb" ? "heatMap" : type) + "]").attr("src");
                var showAttr = $(_this).parents(".mapping_right").attr("showAttr");
                $(_this).parents(".searchResult").find("div .mapping_left").find("li[showAttr=" + showAttr + "]").find("a").css("color", "#2CD869");
                $(_this).parents(".searchResult").find("img[selectAttr=" + (type == "tgb" ? "heatMap" : type) + "]").attr("src", "http://img.biocloud.cn/cloud_beta_1.0/v_1.3/report/report_loading.gif");

                this.postJSONAsync("/report/dge/drawingMaps",
                    {data: JSON.stringify(idList)}, _drawingMapsCallback(filePath, _this, type,str));
            };

            //画图的回调函数
            function _drawingMapsCallback(filePath, _this, type,str) {
                return function (data) {
                    var showAttr = $(_this).parents(".mapping_right").attr("showAttr");
                    if (data == true) {
                        //appUtil.popbox ("画图任务提交成功，请耐心等候几分钟。");

                        $(_this).parents(".searchResult").find("div .mapping_left").find("li[showAttr=" + showAttr + "]").find("a").css("color", "#2CD869");
                        $(_this).parents(".searchResult").find("img[selectAttr=" + (type == "tgb" ? "heatMap" : type) + "]").attr("src", "http://img.biocloud.cn/cloud_beta_1.0/v_1.3/report/report_loading.gif");
                        console.log("任务提交成功：" + new Date().getTime());

                        if(type == "tonglu"){
                            setTimeout(DrawMapping.isComplete(filePath, _this, type), 5000);
                        }else{
                            setTimeout(DrawMapping.isComplete(filePath + ".png", _this, type), 5000);
                        }
                    } else {
                        appUtil.popbox ("画图提交失败，请重试！");
                        $(_this).parents(".searchResult").find("div .mapping_left").find("li[showAttr=" + showAttr + "]").find("a").css("color", "#7a7c7e");
                        $(_this).parents(".searchResult").find("img[selectAttr=" + (type == "tgb" ? "heatMap" : type) + "]").attr("src", str);
                    }
                };
            };

            //判断任务是否完成
            DrawMapping.isComplete = function (filePath, _this, type) {
                this.postJSONAsync("/report/dge/isComplete", {
                    path: filePath,
                    type: type
                }, _isCompleteCallback(filePath, _this, type));
            };

            //查看任务完成之后 回调函数
            function _isCompleteCallback(filePath, _this, type) {
                var showAttr = $(_this).parents(".mapping_right").attr("showAttr");
                return function (data) {
                    if (data.error == "error") {
                        $(_this).parents(".searchResult").find("div .mapping_left").find("li[showAttr=" + showAttr + "]").find("a").css("color", "#FF6347");
                        $(_this).parents(".searchResult").find("img[selectAttr=" + (type == "tgb" ? "heatMap" : type) + "]").parent().html("<img src='http://img.biocloud.cn/cloud_beta_1.0/reg_wrong.png' selectAttr='" + (type == "tgb" ? "heatMap" : type) + "' style='display: inline-block;max-width: none;'/>");
                        appUtil.popbox ("绘图程序报错,请联系客服人员解决问题!");
                        return;
                    } else if (data.error == "true") {
                        $(_this).parents(".searchResult").find("div .mapping_left").find("li[showAttr=" + showAttr + "]").find("a").css("color", "#7A7C7E");
                            $(_this).parents(".searchResult").find("img[selectAttr=" + (type == "tgb" ? "heatMap" : type) + "]").parent().html("<img src='" + "/report/report/getImage?filePath=" + filePath + "' selectAttr='" + (type == "tgb" ? "heatMap" : type) + "' style='display: inline-block;max-width: none;'/>");
                        if (type == "tonglu") {
                            $.post("/report/RNANoRefTrans/getKEGGImage", {
                                dir: filePath
                            },function(imagePath){
                                var  imageName=getBase64Images([imagePath]);
                                $(_this).parents(".searchResult").find("img[selectAttr=" + (type == "tgb" ? "heatMap" : type) + "]").parent().html("<img src='" + imageName + "' selectAttr=" + (type == "tgb" ? "heatMap" : type) + ">");
                                $(_this).parents(".searchResult").find("button[markAttr=download" + type + "]").show();
                                $(_this).parents(".searchResult").find("button[markAttr=download" + type + "]").attr("onclick", "" +
                                "if(confirm('你确认要下载该文件夹？')){" +
                                "window.open('/report/RNANoRefTrans/dwonloadZip?zipPath=" + filePath + "/tonglu.zip&path="+filePath+"','_blank')}");
                            });
                        } else {
                            $(_this).parents(".searchResult").find("button[markAttr=download" + (type == "tgb" ? "heatMap" : type) + "]").show();
                            $(_this).parents(".searchResult").find("button[markAttr=download" + (type == "tgb" ? "heatMap" : type) + "]").attr("onclick", "" +
                            "if(confirm('你确认要下载该文件？')){" +
                            "window.location.href='" + "/util/file/downLoadFile?filepath=" + filePath + "&filename=" + type + ".png'}");
                        }
                    } else {
                        console.log("任务没有完成：" + new Date().getTime());
                        setTimeout(function () {
                            DrawMapping.isComplete(filePath, _this, type);
                        }, 5000);
                    }
                };
            };

            return DrawMapping;
        })();

        var PageView = (function () {
            var PageView = new Object();
            _extends(PageView, Util);
            //切换tab
            PageView.Switch = function (_this) {
                $(_this).addClass("active");
                $(_this).siblings().removeClass("active");
                var gUID = $(_this).attr("Relation");
                $(".rna_no_ref_report_tab").find("div[Relation=" + gUID + "]").show();
                $(".rna_no_ref_report_tab").find("div[Relation=" + gUID + "]").siblings().hide();
            };
            //移除tab
            PageView.removeThisTab = function (_this, table) {
                var gUID = $(_this).parent().parent().attr("Relation");
                //如果移除的是当前显示的tab，则显示其后一个tab，如果后一个tab不存在，则再显示前一个tab。如果都没有则无操作
                if ($(_this).parent().parent().attr("class") == 'active') {
                    var a = $(_this).parent().parent().next().attr("Relation");
                    if (a) {
                        this.Switch($(_this).parent().parent().next());
                    } else if ($(_this).parent().parent().prev().attr("Relation")) {
                        this.Switch($(_this).parent().parent().prev());
                    }
                }
                $("#" + table + " .rna_no_ref_report_tab").find("div[Relation=" + gUID + "]").remove();
                $("#" + table + " " + ".table_data").children("div[class='tab-pane active']").removeClass("active");
                $(_this).parent().parent().remove();
            }

            return PageView;
        })();
        //基因趋势图单选方法---MR
        DataMining.FPKMRadioChange = function (_this) {
            var radio= $(".FPKMValue").find("input[name='Clustering']");
            for (var i=0;i<radio.length;i++){
                if (radio[i].checked){
                    $(radio[i]).parent("div.radio").siblings("ul.row").find("li").find("input[type='text']").removeAttr("disabled");//让他下一个不禁用
                    $(radio[i]).parent("div.radio").siblings("ul.row").find("li").find("input[type='text']").css("background-color", "#ffffff");

                    $(radio[i]).parent("div.radio").parent(".col-sm-12").parent("div.row").siblings("div.row").find("input[type='text']").attr("disabled", "disabled");
                    $(radio[i]).parent("div.radio").parent(".col-sm-12").parent("div.row").siblings("div.row").find("input[type='text']").css("background-color", "#E0E0E0");
                    $(radio[i]).parent("div.radio").parent(".col-sm-12").parent("div.row").siblings("div.row").find("input[type='text']").val("");
                }
            }
        }
        //基因趋势图参数提交
        DataMining.FPKMsubmit = function (_this) {
            if (!appUtil.getPerm($(_this).attr("permDataId"))) {
                appUtil.popbox ("您没有使用的权限,请提升用户等级!");
                return;
            }
            var oriData;
            if(oneOrTwo=="single"){
                oriData=Process.OriginalData.originalDataSin;
            }else if(oneOrTwo=="multiple"){
                oriData=Process.OriginalData.originalDataMul;
            }
            if(countTab("A5E134F315B144B49F0F6343A304DBA0")){
                return ;
            }
            if(!_projectId){
                appUtil.popbox ("没有初始化数据!");
            }

            var a=$("#FPKid").val();
            var b=$("#dajowjdsi").val();
            var c=$("#hldifgusrg").val();
            var num=0;
            if(a==""&&(b==""&&c=="")){
                num=1;
                appUtil.popbox ("输入框不能为空!");
            }
            if(num==0){
                var Kvalue, hierarchicalClustering2, hierarchicalClustering10, key, type;
                $(".FPKMValue li").find("input[type='text']").each(function () {
                    if (this.id == 'dajowjdsi') {
                        hierarchicalClustering2 = this.value.trim();
                    } else if (this.name == 'KvalueClustering') {
                        Kvalue = this.value.trim();
                    } else if (this.id == 'hldifgusrg') {
                        hierarchicalClustering10 = this.value.trim();
                    }
                });
                if (hierarchicalClustering10 && hierarchicalClustering2) {
                    appUtil.popbox ("层次聚类不能输入两种参数");
                    return ;
                }
                var count = this.postJSON("/report/dge/sampleCount",{path:_rootPath}),
                    sample = oriData.size;
                if (Kvalue != null && Kvalue.length > 0) {
                    if (_FPKMvalue2Check(Kvalue.trim(),count,sample)){
                        type = "k-means";
                        key = Kvalue;
                    }
                } else if (hierarchicalClustering2 != null && hierarchicalClustering2.length > 0) {
                    if (_FPKMvalue2Check(hierarchicalClustering2,count,sample)) {
                        type = "K";
                        key = hierarchicalClustering2;
                    }
                } else if (hierarchicalClustering10 != null && hierarchicalClustering10.length > 0) {
                    if (_FPKMvalue10Check(hierarchicalClustering10,count,sample)) {
                        type = "P";
                        key = hierarchicalClustering10;
                    }
                }

                if (key) {
                    $("#FPKMsubmit").hide();
                    $("#FPKMsubmitSpan").append("<img src='http://img.biocloud.cn/AppImg/3.1/025.gif'/>");
                    appUtil.popbox ("基因共表达趋势图绘制中,请您耐心等待!");
                    var json = appUtil.postJSON("/report/dge/FPKMsubmit", {
                        key: key,
                        projectId: _projectId,
                        type: type
                    });
                    _callTrendChart(json,key,type);
                }
            }

        }
        //趋势图回调函数
        function _callTrendChart(json,key,type) {
            var typeName="";
            if(type=="k-means"){
                typeName="K均值聚类-K值:";
            }else if(type=="K"){
                typeName="层次聚类--子树数目:";
            }else{
                typeName="层次聚类--切树比例:";
            }
            var gUID = _addTabToElement("", "基因趋势图", 3, true, false, false, 'A5E134F315B144B49F0F6343A304DBA0');
            var _this = $("#A5E134F315B144B49F0F6343A304DBA0 " + els.MiningCon).find("div[Relation='" + gUID + "']");
            $("#A5E134F315B144B49F0F6343A304DBA0 div[markAttr='"+gUID+"'] .zipDownFileAll").attr("onclick","DgeReport.DataMining.downZip(this,'"+key+"','"+type+"')"); //添加下载zip文件
            $("#A5E134F315B144B49F0F6343A304DBA0 .searchResult[relation='" + gUID + "']").find(".app_experience_directory #title1").text("差异基因挖掘");
            $("#A5E134F315B144B49F0F6343A304DBA0 .searchResult[relation='" + gUID + "']").find(".app_experience_directory #title2").text("基因共表达趋势分析");
            $("#A5E134F315B144B49F0F6343A304DBA0 .searchResult[relation='" + gUID + "']").find(".app_experience_directory #title3").text(typeName+" "+key);
            $("#A5E134F315B144B49F0F6343A304DBA0 .searchResult[relation='" + gUID + "'] .app_experience_flip").hide();//隐藏分页
            $(_this).find(".mining_mapping_cog .img a").html("<img src='http://img.biocloud.cn/cloud_beta_1.0/v_1.3/report/report_loading.gif'>");
            _showContentStyle('A5E134F315B144B49F0F6343A304DBA0');
            setTimeout(DataMining.isTrendSuccess(json, _this, gUID, 'A5E134F315B144B49F0F6343A304DBA0'), 5000);
            appUtil.popbox ("提交后台运行绘图程序,请耐心等待几分钟!");
        }

        //判断趋势图是否完成
        DataMining.isTrendSuccess = function (json, _this, gUID, htmluuid) {
            if (json) {
                appUtil.postJSONAsync("/report/dge/isTrendSuccess", {path: json.jobPath}, _callSuccess(json, _this, gUID, htmluuid));
            }
        }
        //使用java下载zip文件,
        DataMining.downZip = function(_this,key,type){
            if(!($("#FPKMsubmit").is(":hidden"))){
                window.open("/report/RNANoRefTrans/dwonZipFromDir?key="+key+"&type="+type+"&projectId="+_projectId,"_blank");
            }else{
                appUtil.popbox ("正在请求数据,不能下载,请耐心等待几分钟!");
            }

        }
        //判断成功后回调函数
        function _callSuccess(json, _this, gUID, htmluuid) {
            return function (data) {
                if (data.isSuccess == true) {
                    $("#warn_tip_div").remove();
                    var strHtml = "";
                    for (var i = 0; i < data.png.length; i++) {
                        strHtml += ('<li ' + (i == 0 ? 'class="one" style=";background:rgb(219,214,106);"' : "") + ' >' +
                        '<a href="javascript:void(0);" name="' + data.png[i].trim() + '" onclick="DgeReport.DataMining.changehierarchy(this,\'' + gUID + '\')" >' +
                        '<img style=" height: 80px;width: 100px;margin-top: 5px;" src="/report/report/showjobimgByPath?filepath=' + data.png[i].trim() + '"></a></li>');

                    }
                    $(_this).find(".mining_mapping_cog .min_img").html(strHtml);
                    _dataToHTML(data.png[0].trim(), gUID, htmluuid);
                    $(_this).find(".mining_mapping_cog .img_data .m-t-lg a").html("<img style='display: inline-block;' src='" + _projectName + "/report/report/showjobimgByPath?filepath=" + data.png[0].trim() + "'>");
                    $("#A5E134F315B144B49F0F6343A304DBA0 .searchResult[relation='" + gUID + "'] .app_experience_flip").show();//隐藏分页
                    $("#FPKMsubmitSpan").find("img").hide();
                    $("#FPKMsubmit").show();
                } else {
                    console.log("任务没有完成：" + new Date().getTime());
                    setTimeout(function () {
                        DataMining.isTrendSuccess(json, _this, gUID, htmluuid);
                    }, 5000);
                }
            }
        }

        function _dataToHTML(path, gUID, htmluuid) {
            var datafile = path.replace(".png", "").trim();
            var json = appUtil.postJSON("/report/dge/fromNameToData", {path: datafile});
            json.data = _formatDegData(json.data);
            $("#" + htmluuid + " .searchResult[relation='" + gUID + "']").find(".pageBtn").attr("fileName", json.fileName);
            $("#" + htmluuid + " .searchResult[relation='" + gUID + "']").find("div.js-dataMining-ul-download").attr("fileName", json.fileName);
            $("#" + htmluuid + " .searchResult[relation='" + gUID + "']").find(".pageContent").text(json.page);
            $("#" + htmluuid + " .searchResult[relation='" + gUID + "']").find(".pagenow_text").text("1");
            $('#' + htmluuid + ' #data' + gUID).html("");
            $(els.table).tmpl(json).appendTo('#' + htmluuid + ' #data' + gUID);
        }

        //选择趋势图效果
        DataMining.changehierarchy = function (_this, gUID) {
            var parents = $("#A5E134F315B144B49F0F6343A304DBA0 " + els.MiningCon).find("div[Relation='" + gUID + "']");
            $(_this).parent("li").siblings("li").each(function () {
                $(this).css("background", "#ffffff");
            });
            var img = $(_this).attr("name");
            _dataToHTML(img, gUID, "A5E134F315B144B49F0F6343A304DBA0");
            $(_this).parent("li").css("background", "#DBD66A");
            $(_this).find(".mining_mapping_cog .download a").attr("href", "/util/file/downLoadFile?filepath=" + img + "&filename=" + img.split("/")[img.split("/").length - 1]);
            $(parents).find(".mining_mapping_cog .img_data a").html("<img src='/report/report/showjobimgByPath?filepath=" + img + "'>");
        }

        //当层叠聚类输入框数去焦点时input
        DataMining.hierarchyblur = function (_this) {
            if (_this.value) {
                //$(this).parent().append("<span><input name=\"FDR\" type=\"text\"/></span>");
                $(_this).parent().siblings().find("input[type='text']").css("background-color", "#ffffff");
                $(_this).css("background-color", "#ffffff");
                $(_this).parent().siblings().find("input[type='text']").val("");
            } else {
                $(_this).parent().siblings().find("input[type='text']").css("background-color", "#ffffff");
                $(_this).parent().siblings().find("input[type='text']").removeAttr("disabled");
                $(_this).parent().siblings().find("input[type='text']").val("");
            }
            if (_this.value) {
                //$(this).parent().append("<span><input name=\"FDR\" type=\"text\"/></span>");
                $(_this).parent(".li_radioDiv").siblings(".li_radioDiv").find("input[type='text']").css("background-color", "#ffffff");
                $(_this).css("background-color", "#ffffff");
                $(_this).parent(".li_radioDiv").siblings(".li_radioDiv").find("input[type='text']").val("");
            } else {
                $(_this).parent(".li_radioDiv").siblings(".li_radioDiv").find("input[type='text']").css("background-color", "#ffffff");
                $(_this).parent(".li_radioDiv").siblings(".li_radioDiv").find("input[type='text']").removeAttr("disabled");
                $(_this).parent(".li_radioDiv").siblings(".li_radioDiv").find("input[type='text']").val("");
            }
        }
        //基因趋势图输入值校验
        function _FPKMvalue2Check(string,count,sample) {
            string = string.trim();
            var reg = new RegExp('^[0-9]+$','g');
            if (reg.test(string)) {
                if (parseInt(string) > 1 && parseInt(string)< parseInt(count.count.toString()) && parseInt(string) < (Math.pow(3,parseInt(sample)))) {
                    return string;
                }else{
                    appUtil.popbox ("必须输入2以上的整数,最大不超过3的n次方（n=样品个数或者基因个数,即"+(Math.pow(3,parseInt(sample))));
                    return false;
                }
            }else{
                appUtil.popbox ("只能输入整数");
                return false;
            }
        }

        function _FPKMvalue10Check(string,count,sample) {
            string = string.trim();
            var reg = new RegExp('^[0-9]+$','g');
            if (reg.test(string)) {
                if (parseInt(string) > 0 && string < 100 && parseInt(string)<parseInt(count.count.toString()) && parseInt(string) < (Math.pow(3,parseInt(sample)))) {
                    return string;
                }else{
                    appUtil.popbox ("必须输入0-100之间的整数.");
                    return false;
                }
            }else{
                appUtil.popbox ("只能输入数字");
                return false;
            }
        }

        //功能注释搜索回调函数
        function _annoCallback(keywords) {
            return function (data) {
                if(data==""){
                    appUtil.popbox ("没有搜索到结果，请重试！");
                }else{
                    if (data && data[0].data.data.length > 1) {
                        for (var i = 0, l = data.length; i < l; i++) {
                            var json = data[i].data;
                            json.data = _formatData(json.data);
                            var gUID = _addTabToElement("", "关键字搜索", 1, false, true, false, 'A0253C22EF6447C6ACE32D93D81528BD');
                            $('#data' + gUID).html("");
                            $("#A0253C22EF6447C6ACE32D93D81528BD .searchResult[relation='" + gUID + "']").find(".pageBtn").attr("fileName", data[0].data.fileName);
                            $("#A0253C22EF6447C6ACE32D93D81528BD .searchResult[relation='" + gUID + "']").attr("fileName", data[0].data.fileName);
                            $("#A0253C22EF6447C6ACE32D93D81528BD .searchResult[relation='" + gUID + "']").find("div.js-dataMining-ul-download").attr("filename", data[0].data.fileName);
                            $("#A0253C22EF6447C6ACE32D93D81528BD .searchResult[relation='" + gUID + "']").find(".js-dataMining-list-download").attr("fileName", data[0].data.fileName);
                            $("#A0253C22EF6447C6ACE32D93D81528BD .searchResult[relation='" + gUID + "']").find(".js-dataMining-sequence-download").attr("fileName", data[0].data.fileName);
                            $("#A0253C22EF6447C6ACE32D93D81528BD .searchResult[relation='" + gUID + "']").find(".pageContent").text(data[0].data.page);
                            $("#A0253C22EF6447C6ACE32D93D81528BD .searchResult[relation='" + gUID + "']").attr("type", "annoKey");
                            $("#A0253C22EF6447C6ACE32D93D81528BD .searchResult[relation='" + gUID + "'] .app_experience_flip").hide();//隐藏分页
                            $("#A0253C22EF6447C6ACE32D93D81528BD .searchResult[relation='" + gUID + "']").attr("name", keywords);
                            $("#A0253C22EF6447C6ACE32D93D81528BD .searchResult[relation='" + gUID + "']").find(".app_experience_directory #title1").text("基因挖掘");
                            $("#A0253C22EF6447C6ACE32D93D81528BD .searchResult[relation='" + gUID + "']").find(".app_experience_directory #title2").text("基因信息检索");
                            $("#A0253C22EF6447C6ACE32D93D81528BD .searchResult[relation='" + gUID + "']").find(".app_experience_directory #title3").text("类型:功能关键词,参数:"+keywords);
                            $(els.table).tmpl(json).appendTo('#A0253C22EF6447C6ACE32D93D81528BD #data' + gUID);
                        }
                        _showContentStyle('A0253C22EF6447C6ACE32D93D81528BD');
                    } else {
                        appUtil.popbox ("没有搜索到结果，请重试！");
                    }
                }
                $("#annoSearchKeywordSpan").find("img").hide();
                $("#annoSearchKeyword").show();
            };
        };

        //功能id序列搜索回调函数
        function _annoCallbackId(keywords) {
            return function (data) {
                if(JSON.stringify(data)=="{}"){
                    appUtil.popbox ("没有搜索到结果，请重试！");
                }else{
                    data.data = _formatData(data.data);
                    var changeValue = $(els.searchType).val();
                    if (changeValue == "SequenceSearch") {
                        var gUID = _addTabToElement("", "序列搜索", 1, false, true, false, 'A0253C22EF6447C6ACE32D93D81528BD',false,"default");
                        //var keywords = data.searchName;
                    } else {
                        var gUID = _addTabToElement("", "基因ID搜索", 1, false, true, false, 'A0253C22EF6447C6ACE32D93D81528BD',false,"default");
                    }
                    $('#A0253C22EF6447C6ACE32D93D81528BD #data' + gUID).html("");
                    $(els.table).tmpl(data).appendTo('#A0253C22EF6447C6ACE32D93D81528BD #data' + gUID);
                    $("#A0253C22EF6447C6ACE32D93D81528BD .searchResult[relation='" + gUID + "']").find(".pageBtn").attr("fileName", data.fileName);
                    $("#A0253C22EF6447C6ACE32D93D81528BD .searchResult[relation='" + gUID + "']").find("div.js-dataMining-ul-download").attr("filename", data.fileName);
                    $("#A0253C22EF6447C6ACE32D93D81528BD .searchResult[relation='" + gUID + "']").find(".pageContent").text(data.page);
                    $("#A0253C22EF6447C6ACE32D93D81528BD .searchResult[relation='" + gUID + "']").attr("type", "SequenceSearch");
                    $("#A0253C22EF6447C6ACE32D93D81528BD .searchResult[relation='" + gUID + "']").attr("name", keywords);
                    $("#A0253C22EF6447C6ACE32D93D81528BD .searchResult[relation='" + gUID + "']").find(".app_experience_directory #title1").text("基因挖掘");
                    $("#A0253C22EF6447C6ACE32D93D81528BD .searchResult[relation='" + gUID + "']").find(".app_experience_directory #title2").text("基因信息检索");
                    $("#A0253C22EF6447C6ACE32D93D81528BD .searchResult[relation='" + gUID + "'] .app_experience_flip").hide();//隐藏分页
                    if(changeValue == "SequenceSearch"){
                        $("#A0253C22EF6447C6ACE32D93D81528BD .searchResult[relation='" + gUID + "']").find(".app_experience_directory #title3").text("类型:基因序列片段,参数:"+keywords);
                    }else{
                        $("#A0253C22EF6447C6ACE32D93D81528BD .searchResult[relation='" + gUID + "']").find(".app_experience_directory #title3").text("类型:基因编号列表,参数:"+keywords);
                    }
                    $("#A0253C22EF6447C6ACE32D93D81528BD .searchResult[relation='" + gUID + "']").attr("name", keywords);
                    _showContentStyle('A0253C22EF6447C6ACE32D93D81528BD');
                }
                $("#annoSearchKeywordSpan").find("img").hide();
                $("#annoSearchKeyword").show();
            };
        };

        //输入框为空校验
        function _checkNull(value){
            if(value.trim().length==0){
                appUtil.popbox ("输入框不能为空！");
                throw "";
            }
        }

        //注释功能和序列搜索、差异结果查询、差异表达基因Venn、基因共表达趋势分析、样品基因表达Venn:分页模块
        var PageContent = (function () {
            var PageContent = new Object();
            _extends(PageContent, Util);
            //首页
            PageContent.first = function (_this) {
                if (!checkPage("first", _this)) {
                    appUtil.popbox ("已经到达首页!");
                    return;
                }
                flipPage(1, _this);
            }
            //上一页
            PageContent.upPage = function (_this) {
                if (!checkPage("upPage", _this)) {
                    appUtil.popbox ("已经到达首页!");
                    return;
                }
                var nowPage = $(_this).parents().parents("ul").find(".pagenow_text").text();
                flipPage(parseInt(nowPage) - 1, _this);
            }
            //下一页
            PageContent.downPage = function (_this) {
                if (!checkPage("downPage", _this)) {
                    appUtil.popbox ("已经到达尾页!");
                    return;
                }
                var nowPage = $(_this).parents().parents("ul").find(".pagenow_text").text();
                flipPage(parseInt(nowPage) + 1, _this);
            }
            //尾页
            PageContent.end = function (_this) {
                if (!checkPage("end", _this)) {
                    appUtil.popbox ("已经到达尾页!");
                    return;
                }
                var maxPage = $(_this).parents().parents("ul").find(".pageContent").text();
                flipPage(parseInt(maxPage), _this);
            }
            //手动跳转
            PageContent.change = function (_this) {
                if (!checkPage("change", _this)) {
                    appUtil.popbox ("选择超出范围");
                    $(_this).val("1");
                    return;
                }
                var nowPage = $(_this).parents().parents("ul").find(".pagenow_text").text();
                flipPage(parseInt(nowPage), _this);
            }
            //提交翻页
            function flipPage(page, _this) {
                var gUID = $(_this).parents(".pageBtn").attr("gUID"),
                    uuid = $(_this).parents(".pageBtn").attr("fileName"),
                    htmluuid = $(_this).parents(".mining").attr("id");
                var data = appUtil.postJSON("/report/dge/dataContentPage", {
                    page: page,
                    uuid: uuid
                });
                if (data != null) {
                    $(_this).parents("ul").find(".pagenow_text").text(page);
                    data.data = _formatDegData(data.data);
                    $('#' + htmluuid + ' #data' + gUID).html("");
                    $(els.table).tmpl(data).appendTo('#' + htmluuid + ' #data' + gUID);
                }
            }

            //判断是否可以翻页
            function checkPage(direction, _this) {
                var nowPage = $(_this).parents().parents("ul").find(".pagenow_text").text(),
                    maxPage = $(_this).parents().parents("ul").find(".pageContent").text();

                if (parseInt(nowPage) < 0 && parseInt(nowPage) > parseInt(maxPage)) {
                    return false;
                }
                switch (direction) {
                    case "upPage":
                        if (nowPage != '1' && nowPage != null) {
                            return true;
                        }
                        break;
                    case "first":
                        if (nowPage != '1' && nowPage != null) {
                            return true;
                        }
                        break;
                    case "downPage":
                        if (nowPage != maxPage && nowPage != null) {
                            return true;
                        }
                        break;
                    case "end":
                        if (nowPage != maxPage && nowPage != null) {
                            return true;
                        }
                        break;
                    case "change":
                        var reg = new RegExp(/^[0-9]+$/g);
                        if (!reg.test(nowPage)) {
                            appUtil.popbox ("请输入正确的数字");
                            return;
                        }
                        if (parseInt(nowPage) > 0 && parseInt(nowPage) <= parseInt(maxPage) && nowPage != null) {
                            return true;
                        }
                }
                return false;
            }

            return PageContent;
        })();

        //蛋白互作图提交,
        DataMining.proteinInteraction = function (_this) {
            var oriData;
            if(oneOrTwo==""){
                if(typeselect=="single"){
                    oriData=Process.OriginalData.originalDataSin;
                }else if(typeselect=="multiple"){
                    oriData=Process.OriginalData.originalDataMul;
                }
            }else if(oneOrTwo=="single"){
                oriData=Process.OriginalData.originalDataSin;
            }else if(oneOrTwo=="multiple"){
                oriData=Process.OriginalData.originalDataMul;
            }
            var fileName = $(_this).parents("div.js-dataMining-ul-download").attr("fileName");
            if (isCytoscape) {
                _proteinInteraction(fileName);
                return;
            }
            var json = appUtil.postJSON("/report/dge/isCytoscape", {projectId: _projectId});
            if (json != null) {
                if (json.error == "true") {
                    isCytoscape = true;
                    _proteinInteraction(fileName);
                    return;
                } else {
                    if (oriData.size == 1) {
                        appUtil.popbox ("项目紧包含一个样品无法制作互作图");
                    }
                    var data = appUtil.postJSON("/report/dge/cytoscape", {projectId: _projectId});
                    if (data != null) {
                        if (data.error == 'true')
                            appUtil.popbox ("提交蛋白图制作大概需要1到6个小时运行程序,可以在'我的任务'中查看进度(占用任务数资源),请耐心等待!");
                        else
                            appUtil.popbox ("任务提交失败:" + data.message + "");
                    }
                }
            }
        }
        //蛋白互作图提交程序
        function _proteinInteraction(fileName) {
            window.open("/report/dge/proteinInteraction?projectId=" + _projectId + "&fileName=" + fileName, "_blank");
        }

        DataMining.PageContent = PageContent;
        DataMining.PageView = PageView;
        DataMining.DegGroupSelect = DegGroupSelect;
        DataMining.DrawMapping = DrawMapping;
        return DataMining;
    })();
    //蛋白互作图
    ProteinIntMapping = (function () {
        var ProteinIntMapping = new Object();
        _extends(ProteinIntMapping, Util);

        //加载项目差异分析
        ProteinIntMapping.init = function () {
            var SepDefined,ComDefined;
            if(("Sep" in Process.detailCfgNext)){
                SepDefined = Process.detailCfgNext.Sep;
            }
            if(!("Com" in Process.detailCfgNext)){
                for(var item in Process.detailCfgNext){
                    if(item == "Com"){
                        ComDefined += (Process.detailCfgNext.item);
                    }
                }
                if(ComDefined!=null ){
                    ComDefined = ComDefined.substring(0,ComDefined.length-1);
                }
            }
            if (ComDefined == null && SepDefined == null) {
                $("#dnahfishofuahwefnskf").html("<option value='null'>(NULL)</option>");
            } else {
                var strHtml = "<option value=''>请选择差异分析</option>";
                if (ComDefined != null) {
                    var com = ComDefined.split("!");
                    for (var c in com) {
                        strHtml += "<option value='" + com[c].replace(",", "_vs_") + "'>" + com[c].replace(",", "_vs_") + "</option>";
                    }
                }
                if (SepDefined != null) {
                    var sep = SepDefined.split("!");
                    for (var s in sep) {
                        var sre = sep[s].replace(/,/g, "_");
                        sre = sre.replace(";", "_vs_");
                        strHtml += "<option value='" + sre + "'>" + sre + "</option>";
                    }
                }
                $("#dnahfishofuahwefnskf").html(strHtml);
            }
        }
        //蛋白图后台请求
        ProteinIntMapping.differenceAnalysis = function (_this) {
            var oriData;
            if(oneOrTwo==""){
                if(typeselect=="single"){
                    oriData=Process.OriginalData.originalDataSin;
                }else if(typeselect=="multiple"){
                    oriData=Process.OriginalData.originalDataMul;
                }
            }else if(oneOrTwo=="single"){
                oriData=Process.OriginalData.originalDataSin;
            }else if(oneOrTwo=="multiple"){
                oriData=Process.OriginalData.originalDataMul;
            }
            if (isCytoscape == false) {
                if (oriData.size == 1) {
                    appUtil.popbox ("项目紧包含一个样品无法制作互作图");
                }
                var data = appUtil.postJSON(_projectName + "/refReport/RNARefReport/cytoscape", {projectId: _projectId});
                if (data != null) {
                    if (data.error == 'true')
                        appUtil.popbox ("提交蛋白图制作大概需要1到6个小时运行程序,可以在'我的任务'中查看进度(占用任务数资源),请耐心等待!");
                    else
                        appUtil.popbox ("任务提交失败:" + data.message);
                }
            } else {
                var change = $("#dnahfishofuahwefnskf").val();
                if (change == null || change == "") {
                    appUtil.popbox ("请选择差异分组!");
                }
                var data = appUtil.postJSON(_projectName + "/refReport/RNARefReport/proteinIntMapping", {
                    projectId: _projectId,
                    change: change
                });
                if (data != null) {
                    if (data.nodesdata) {
                        $("#aodiajwodajoiwdj").attr("name", data.change);
                        _callProteinIntList(data);
                        _callProteinInt(data);
                    }
                } else {
                    appUtil.popbox ("没有相对应的蛋白图文件");
                }
            }
        };
        //蛋白互作图回调函数去画图
        function _callProteinInt(data) {
            if (!data) {
            } else {
                var nodesdata = data.nodesdata;
                var edgesdata = data.edgesdata;
                var layout = $("#gjdfjgidhusgf").val();
                cytoscape({
                    container: $('#daowjdoajiwdojifhdurg')[0],
                    style: cytoscape.stylesheet()
                        //节点
                        .selector('node')
                        .css({
                            'content': 'data(id)',
                            'text-valign': 'center',//节点形状
                            'color': '#000000',//文字颜色
                            'background-color': '#DA9494',//背景色
                            'width': 'mapData(score, 0, 100, 10, 200)',
                            'height': 'mapData(score, 0, 100, 10, 200)'
                        })
                        //线
                        .selector('edge')
                        .css({
                            'line-color': '#114ED0',
                            'width': 2
                        })
                        //被选中的节点
                        .selector(':selected')
                        .css({
                            'background-color': '#0D882F',//被选择节点颜色
                            'line-color': '#777777'//被选中线的颜色
                        }),
                    elements: {
                        nodes: nodesdata,
                        edges: edgesdata
                    },
                    layout: {
                        name: layout == "" ? 'concentric' : layout,//'concentric',
                        padding: 5,
                        radius: 100,
                        levelWidth: function (nodes) {
                            return 0.1;
                        }
                    }
                });
                $("#mirhugsidngkixncvbuihsf").show();
            }
        }

        //蛋白互作图中参数集合
        function _callProteinIntList(data) {
            if (data) {
                var edgesdata = data.edgesdata, strHtml = "";
                for (var i in edgesdata) {
                    var edge = edgesdata[i];
                    strHtml += "<p>" + edge.data.source + "--" + edge.data.target + "</p>";
                }
                $("#daojfosdfosndfoef").attr("name", data.change);

                $("#djawoihfiushfse").html(strHtml);
            }
        }

        //下载蛋白图数据
        ProteinIntMapping.downloadProteinInt = function (_this) {
            var name = $("#daojfosdfosndfoef").attr("name");
            if (name) {
                window.location.href = _projectName + '/refReport/RNARefReport/downloadData?name=' + name + '&projectId=' + _projectId;
            }
        }
        //另起一页放大图片
        ProteinIntMapping.changeResolution = function (_this) {
            var change = _this.name;
            window.open(_projectName + "/refReport/RNARefReport/networkdiagram?projectId=" + _projectId + "&change=" + change, "_blank");
        }
        //下载图片
        ProteinIntMapping.convertCanvasToImage = function (_this) {
            var canvas = null;
            $("#daowjdoajiwdojifhdurg div").find("canvas").each(function () {
                if ($(this).css("z-index") == 1) {
                    canvas = this;
                }
            });
            var image = new Image();
            var canvas2 = renderScreenshot(canvas, "4");
            image.src = canvas2;
            image.src = image.src.replace("image/png", 'image/octet-stream');
            var saveFile = function (data, filename) {
                var save_link = document.createElementNS('http://www.w3.org/1999/xhtml', 'a');
                save_link.href = data;
                save_link.download = filename;

                var event = document.createEvent('MouseEvents');
                event.initMouseEvent('click', true, false, window, 0, 0, 0, 0, 0, false, false, false, false, 0, null);
                save_link.dispatchEvent(event);
            };
            var filename = 'appUtil_' + (new Date()).getTime() + '.' + "png";
            saveFile(image.src, filename);
        }
        function renderScreenshot(canvas, scaleFactor) {
            var ctx = canvas.getContext('2d');
            var screenshotCanvas = document.createElement('canvas');
            var screenshotCtx = screenshotCanvas.getContext('2d');
            screenshotCanvas.width = canvas.width * scaleFactor;
            screenshotCanvas.height = canvas.height * scaleFactor;
            screenshotCtx.drawImage(canvas, 0, 0, canvas.width * scaleFactor, canvas.height * scaleFactor);
            return screenshotCanvas.toDataURL();
        }

        return ProteinIntMapping;
    })();
    //页面样式模块
    PageStyle = (function () {
        PageStyle = PageStyle || {};
        _extends(PageStyle, Util);
        var initialCode = PageStyle.initialCode = new Object(),
            callCode = PageStyle.callCode = new Object();

        initialCode.topTab = function () {
            //顶层tab切换，流程定制，数据挖掘
            $(".repnav").delegate("li", "click", function () {
                $(".pop_up").remove();
                $(this).addClass("one");
                var help = $(this).attr("name");
                $("." + help).siblings().attr("style", "display:none;");
                $("." + help).show();
                $(this).siblings().removeClass("one");
                var num = $(this).attr("num");
                $(".qweqwe" + num).show();
                var brothers = $(this).siblings();
                for (var i = 0; i < brothers.length; i++) {
                    var numm = $(brothers[i]).attr("num");
                    $(".qweqwe" + numm).hide();
                }
            });
        };
        initialCode.originalDataLiPop = function (_this, hover) {
            if (hover == 'show') {
                $(_this).next().show();
            } else {
                $(_this).next().hide();
            }
        }
        //帮助弹出框
        initialCode.helpPop = function () {
            $(".problem a").hover(function () {
                $(this).siblings(".problem_content").show();
            }, function () {
                $(".problem_content").hover(function(){
                    $(this).show();
                },function(){
                    $(".problem a").siblings(".problem_content").hide();
                });
                $(this).siblings(".problem_content").hide();
            });
        };


        //当用户选择自定义参数的时候，该下拉框修改成输出框
        initialCode.CustomFDR = function () {
            $(".morw select[name=FDR]").change(function () {
                var value = $(this).val();
                if (value == 0) {
                    $(this).parent().append("<span><input name=\"FDR\" type=\"text\"/></span>");
                    $(this).remove();
                }
            });
        };

        //当用户选择自定义参数的时候，该下拉框修改成输出框
        initialCode.CustomFC = function () {
            $(".morw select[name=FC]").change(function () {
                var value = $(this).val();
                if (value == 0) {
                   $(this).parent().append("<span><input name=\"FC\" type=\"text\"/></span>");
                    $(this).remove();
                }
            });
        };

        //只上调和只下调 只能选择一个
        initialCode.degSelect = function () {
            $("input[name=degSelect]").click(function () {
                var isSelect = $(this).is(':checked');
                if (isSelect === true) {
                    $("input[name=degSelect]").prop("checked", false);
                    $(this).prop("checked", true);
                }
            });
        };

        //只上调和只下调 只能选择一个
        initialCode.degGrouping = function () {
            $("input[name=degGroupingRegulated]").click(function () {
                var isSelect = $(this).is(':checked');
                if (isSelect === true) {
                    $("input[name=degGroupingRegulated]").prop("checked",false);
                    $(this).prop("checked", true);
                }
            });
        };

        //数据挖掘内容加载到页面中是，给画图tab添加切换事件
        callCode.zuotu = function () {
            $(".mapping_left").delegate("li", "click", function () {
                $(this).addClass("active");
                $(this).siblings().removeClass("active");
                var qwe = $(this).attr("showAttr");
                $(".mapping_right[showAttr=" + qwe + "]").show();
                var qwes = $(this).siblings();
                for (var i = 0; i < qwes.length; i++) {
                    var otherClass = $(qwes[i]).attr("showAttr");
                    $(".mapping_right[showAttr=" + otherClass + "]").hide();
                }
            })
        };

        //差异分组选择限制
        callCode.degGroCheckbox = function (_this) {
            var degVal = $(_this).val(), name = $(_this).attr("name") == "Contrast_check" ? "Experiment_check" : "Contrast_check";
            if ($(_this).is(":checked")) {
                $(_this).parents().parents().parents().parents(".div_differencesGroup").siblings(".div_differencesGroup").find("li input[name='" + name + "'][value='" + degVal + "']").attr("disabled", "disabled");
            } else {
                $(_this).parents().parents().parents().parents(".div_differencesGroup").siblings(".div_differencesGroup").find("li input[name='" + name + "'][value='" + degVal + "']").removeAttr("disabled");
            }
        }

        //用来移除liDOM
        callCode.removeLi = function (_this) {
            $(_this).parent().parent().remove();
        };

        //用来关闭弹出框
        callCode.closePop = function () {
            $(".modal").remove();
            $(".modal-backdrop").hide();
        };

        //关闭本页
        callCode.custom_close = function () {
            if (confirm("您确定要离开本页吗？")) {
                window.opener = null;
                window.open('', '_self');
                window.close();
            }
        };

        return PageStyle;

    })();
    //挂载几个模块到Dge下面
    Dge.Process=Process;
    Dge.PageStyle=PageStyle;
    Dge.DataMining = DataMining;
    Dge.ProteinIntMapping = ProteinIntMapping;

})(window,$,appUtil);

$(document).ready(function(){

    $("#wizard").steps();
    $("#form").steps({
        bodyTag: "fieldset",
        onStepChanging: function (event, currentIndex, newIndex) {
            // Always allow going backward even if the current step contains invalid fields!
            if (currentIndex > newIndex) {
                return true;
            }

            // Forbid suppressing "Warning" step if the user is to young
            if (newIndex === 3 && Number($("#age").val()) < 18) {
                return false;
            }

            var form = $(this);

            // Clean up if user went backward before
            if (currentIndex < newIndex) {
                // To remove error styles
                $(".body:eq(" + newIndex + ") label.error", form).remove();
                $(".body:eq(" + newIndex + ") .error", form).removeClass("error");
            }

            // Disable validation on fields that are disabled or hidden.
            form.validate().settings.ignore = ":disabled,:hidden";

            // Start validation; Prevent going forward if false
            return form.valid();
        },
        onStepChanged: function (event, currentIndex, priorIndex) {
            // Suppress (skip) "Warning" step if the user is old enough.
            if (currentIndex === 2 && Number($("#age").val()) >= 18) {
                $(this).steps("next");
            }

            // Suppress (skip) "Warning" step if the user is old enough and wants to the previous step.
            if (currentIndex === 2 && priorIndex === 3) {
                $(this).steps("previous");
            }
        },
        onFinishing: function (event, currentIndex) {
            var form = $(this);

            // Disable validation on fields that are disabled.
            // At this point it's recommended to do an overall check (mean ignoring only disabled fields)
            form.validate().settings.ignore = ":disabled";

            // Start validation; Prevent form submission if false
            return form.valid();
        },
        onFinished: function (event, currentIndex) {
            var form = $(this);

            // Submit form input
            form.submit();
        }
    }).validate({
        errorPlacement: function (error, element) {
            element.before(error);
        },
        rules: {
            confirm: {
                equalTo: "#password"
            }
        }
    });
    $('.jstree1').jstree({
        'core': {
            'check_callback': true
        },
        'plugins': ['types', 'dnd'],
        'types': {
            'default': {
                'icon': 'fa fa-folder'
            },
            'html': {
                'icon': 'fa fa-file-code-o'
            },
            'svg': {
                'icon': 'fa fa-file-picture-o'
            },
            'css': {
                'icon': 'fa fa-file-code-o'
            },
            'img': {
                'icon': 'fa fa-file-image-o'
            },
            'js': {
                'icon': 'fa fa-file-text-o'
            }

        }
    });

    var projectId;
    var isBlank = $("#isBlank").val();
    var rootPath = $("#rootPath").val();
    var contextpath_hidden = $("#contextpath_hidden").val();
    var drawMapPath = "/share/bioCloud/cloud/temp/";

    if (isBlank) {
        projectId=$("#projectId").val();
    }/*else{
        var array = window.location.href.split("/");
        projectId = array[array.length - 1].split("#")[0].split("?")[0];
    }*/

    $(".help_top_no_ref").find("#home_help").removeAttr("style");
    DgeReport.init({
        rootPath: rootPath,
        isBlank: isBlank ? false : true,
        projectName: contextpath_hidden,
        drawMapPath: drawMapPath,
        projectId: projectId
    }).run();
});