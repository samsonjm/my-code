<!DOCTYPE html>
<html class="" lang="en">
<head prefix="og: http://ogp.me/ns#">
<meta charset="utf-8">
<meta content="IE=edge" http-equiv="X-UA-Compatible">
<meta content="object" property="og:type">
<meta content="GitLab" property="og:site_name">
<meta content="README.md · master · Samson, Jonathan / runhisat" property="og:title">
<meta content="Runs a HISAT2 pipeline to prepare RNAseq data for analysis using DESeq." property="og:description">
<meta content="https://git.wur.nl/assets/gitlab_logo-7ae504fe4f68fdebb3c2034e36621930cd36ea87924c11ff65dbcb8ed50dca58.png" property="og:image">
<meta content="64" property="og:image:width">
<meta content="64" property="og:image:height">
<meta content="https://git.wur.nl/samso008/runhisat/blob/master/README.md" property="og:url">
<meta content="summary" property="twitter:card">
<meta content="README.md · master · Samson, Jonathan / runhisat" property="twitter:title">
<meta content="Runs a HISAT2 pipeline to prepare RNAseq data for analysis using DESeq." property="twitter:description">
<meta content="https://git.wur.nl/assets/gitlab_logo-7ae504fe4f68fdebb3c2034e36621930cd36ea87924c11ff65dbcb8ed50dca58.png" property="twitter:image">

<title>README.md · master · Samson, Jonathan / runhisat · GitLab</title>
<meta content="Runs a HISAT2 pipeline to prepare RNAseq data for analysis using DESeq." name="description">
<link rel="shortcut icon" type="image/png" href="/assets/favicon-7901bd695fb93edb07975966062049829afb56cf11511236e61bcf425070e36e.png" id="favicon" data-original-href="/assets/favicon-7901bd695fb93edb07975966062049829afb56cf11511236e61bcf425070e36e.png" />
<link rel="stylesheet" media="all" href="/assets/application-def1880ada798c68ee010ba2193f53a2c65a8981871a634ae7e18ccdcd503fa3.css" />
<link rel="stylesheet" media="print" href="/assets/print-74c3df10dad473d66660c828e3aa54ca3bfeac6d8bb708643331403fe7211e60.css" />


<link rel="stylesheet" media="all" href="/assets/highlight/themes/white-3144068cf4f603d290f553b653926358ddcd02493b9728f62417682657fc58c0.css" />
<script>
//<![CDATA[
window.gon={};gon.api_version="v4";gon.default_avatar_url="https://git.wur.nl/assets/no_avatar-849f9c04a3a0d0cea2424ae97b27447dc64a7dbfae83c036c45b403392f0e8ba.png";gon.max_file_size=10;gon.asset_host=null;gon.webpack_public_path="/assets/webpack/";gon.relative_url_root="";gon.shortcuts_path="/help/shortcuts";gon.user_color_scheme="white";gon.gitlab_url="https://git.wur.nl";gon.revision="2417d5becc7";gon.gitlab_logo="/assets/gitlab_logo-7ae504fe4f68fdebb3c2034e36621930cd36ea87924c11ff65dbcb8ed50dca58.png";gon.sprite_icons="/assets/icons-4009ebf96719f129f954d643e65f87152b5e6c4a4917130c4a696beb54af9949.svg";gon.sprite_file_icons="/assets/file_icons-7262fc6897e02f1ceaf8de43dc33afa5e4f9a2067f4f68ef77dcc87946575e9e.svg";gon.emoji_sprites_css_path="/assets/emoji_sprites-289eccffb1183c188b630297431be837765d9ff4aed6130cf738586fb307c170.css";gon.test_env=false;gon.suggested_label_colors={"#0033CC":"UA blue","#428BCA":"Moderate blue","#44AD8E":"Lime green","#A8D695":"Feijoa","#5CB85C":"Slightly desaturated green","#69D100":"Bright green","#004E00":"Very dark lime green","#34495E":"Very dark desaturated blue","#7F8C8D":"Dark grayish cyan","#A295D6":"Slightly desaturated blue","#5843AD":"Dark moderate blue","#8E44AD":"Dark moderate violet","#FFECDB":"Very pale orange","#AD4363":"Dark moderate pink","#D10069":"Strong pink","#CC0033":"Strong red","#FF0000":"Pure red","#D9534F":"Soft red","#D1D100":"Strong yellow","#F0AD4E":"Soft orange","#AD8D43":"Dark moderate orange"};gon.first_day_of_week=0;gon.ee=false;gon.current_user_id=1291;gon.current_username="lohr001";gon.current_user_fullname="Lohr, Julia";gon.current_user_avatar_url="https://secure.gravatar.com/avatar/6e7f9ff27c82118fa1299609f0124e5f?s=80\u0026d=identicon";
//]]>
</script>


<script src="/assets/webpack/runtime.15f5b13a.bundle.js" defer="defer"></script>
<script src="/assets/webpack/main.f8ac0dda.chunk.js" defer="defer"></script>
<script src="/assets/webpack/commons~pages.admin.clusters~pages.admin.clusters.destroy~pages.admin.clusters.edit~pages.admin.clus~17e57b7f.1e0b3713.chunk.js" defer="defer"></script>
<script src="/assets/webpack/commons~pages.groups.milestones.edit~pages.groups.milestones.new~pages.projects.blame.show~pages.pro~d3e579ac.0dc404e7.chunk.js" defer="defer"></script>
<script src="/assets/webpack/pages.projects.blob.show.d9236e3f.chunk.js" defer="defer"></script>
<script>
//<![CDATA[
window.uploads_path = "/samso008/runhisat/uploads";



//]]>
</script>
<meta name="csrf-param" content="authenticity_token" />
<meta name="csrf-token" content="7a+1h/bK5LuZUXFbQiwkR+/o2xqXnSasfHReWTmGhgCnNzdOXSRtOwQewOpW9BDr7kGU2IztsUlrfT7WTKjEMg==" />

<meta content="origin-when-cross-origin" name="referrer">
<meta content="width=device-width, initial-scale=1, maximum-scale=1" name="viewport">
<meta content="#474D57" name="theme-color">
<link rel="apple-touch-icon" type="image/x-icon" href="/assets/touch-icon-iphone-5a9cee0e8a51212e70b90c87c12f382c428870c0ff67d1eb034d884b78d2dae7.png" />
<link rel="apple-touch-icon" type="image/x-icon" href="/assets/touch-icon-ipad-a6eec6aeb9da138e507593b464fdac213047e49d3093fc30e90d9a995df83ba3.png" sizes="76x76" />
<link rel="apple-touch-icon" type="image/x-icon" href="/assets/touch-icon-iphone-retina-72e2aadf86513a56e050e7f0f2355deaa19cc17ed97bbe5147847f2748e5a3e3.png" sizes="120x120" />
<link rel="apple-touch-icon" type="image/x-icon" href="/assets/touch-icon-ipad-retina-8ebe416f5313483d9c1bc772b5bbe03ecad52a54eba443e5215a22caed2a16a2.png" sizes="152x152" />
<link color="rgb(226, 67, 41)" href="/assets/logo-d36b5212042cebc89b96df4bf6ac24e43db316143e89926c0db839ff694d2de4.svg" rel="mask-icon">
<meta content="/assets/msapplication-tile-1196ec67452f618d39cdd85e2e3a542f76574c071051ae7effbfde01710eb17d.png" name="msapplication-TileImage">
<meta content="#30353E" name="msapplication-TileColor">




</head>

<body class="ui-light  gl-browser-firefox gl-platform-windows" data-find-file="/samso008/runhisat/find_file/master" data-group="" data-page="projects:blob:show" data-project="runhisat">

<script>
//<![CDATA[
gl = window.gl || {};
gl.client = {"isFirefox":true,"isWindows":true};


//]]>
</script>


<header class="navbar navbar-gitlab navbar-expand-sm js-navbar" data-qa-selector="navbar">
<a class="sr-only gl-accessibility" href="#content-body" tabindex="1">Skip to content</a>
<div class="container-fluid">
<div class="header-content">
<div class="title-container">
<h1 class="title">
<a title="Dashboard" id="logo" href="/"><img class="brand-header-logo lazy" data-src="/uploads/-/system/appearance/header_logo/1/wur_logo.png" src="data:image/gif;base64,R0lGODlhAQABAAAAACH5BAEKAAEALAAAAAABAAEAAAICTAEAOw==" />
</a></h1>
<ul class="list-unstyled navbar-sub-nav">
<li id="nav-projects-dropdown" class="home dropdown header-projects qa-projects-dropdown" data-track-label="projects_dropdown" data-track-event="click_dropdown" data-track-value=""><button class="btn" data-toggle="dropdown" type="button">
Projects
<svg class="caret-down"><use xlink:href="/assets/icons-4009ebf96719f129f954d643e65f87152b5e6c4a4917130c4a696beb54af9949.svg#angle-down"></use></svg>
</button>
<div class="dropdown-menu frequent-items-dropdown-menu">
<div class="frequent-items-dropdown-container">
<div class="frequent-items-dropdown-sidebar qa-projects-dropdown-sidebar">
<ul>
<li class=""><a class="qa-your-projects-link" href="/dashboard/projects">Your projects
</a></li><li class=""><a href="/dashboard/projects/starred">Starred projects
</a></li><li class=""><a href="/explore">Explore projects
</a></li></ul>
</div>
<div class="frequent-items-dropdown-content">
<div data-project-id="5460" data-project-name="runhisat" data-project-namespace="Samson, Jonathan / runhisat" data-project-web-url="/samso008/runhisat" data-user-name="lohr001" id="js-projects-dropdown"></div>
</div>
</div>

</div>
</li><li id="nav-groups-dropdown" class="home dropdown header-groups qa-groups-dropdown" data-track-label="groups_dropdown" data-track-event="click_dropdown" data-track-value=""><button class="btn" data-toggle="dropdown" type="button">
Groups
<svg class="caret-down"><use xlink:href="/assets/icons-4009ebf96719f129f954d643e65f87152b5e6c4a4917130c4a696beb54af9949.svg#angle-down"></use></svg>
</button>
<div class="dropdown-menu frequent-items-dropdown-menu">
<div class="frequent-items-dropdown-container">
<div class="frequent-items-dropdown-sidebar qa-groups-dropdown-sidebar">
<ul>
<li class=""><a class="qa-your-groups-link" href="/dashboard/groups">Your groups
</a></li><li class=""><a href="/explore/groups">Explore groups
</a></li></ul>
</div>
<div class="frequent-items-dropdown-content">
<div data-user-name="lohr001" id="js-groups-dropdown"></div>
</div>
</div>

</div>
</li><li class="d-none d-xl-block"><a class="dashboard-shortcuts-activity" href="/dashboard/activity">Activity
</a></li><li class="d-none d-xl-block"><a class="dashboard-shortcuts-milestones" href="/dashboard/milestones">Milestones
</a></li><li class="d-none d-xl-block"><a class="dashboard-shortcuts-snippets qa-snippets-link" href="/dashboard/snippets">Snippets
</a></li>
<li class="d-xl-none dropdown header-more">
<a data-toggle="dropdown" href="#">
More
<svg class="caret-down"><use xlink:href="/assets/icons-4009ebf96719f129f954d643e65f87152b5e6c4a4917130c4a696beb54af9949.svg#angle-down"></use></svg>
</a>
<div class="dropdown-menu">
<ul>
<li class=""><a href="/dashboard/activity">Activity
</a></li><li class=""><a class="dashboard-shortcuts-milestones" href="/dashboard/milestones">Milestones
</a></li><li class=""><a class="dashboard-shortcuts-snippets" href="/dashboard/snippets">Snippets
</a></li>
<li class="dropdown d-lg-none">

</li>
<li class="d-lg-none"><a href="/-/instance_statistics">Instance Statistics
</a></li></ul>
</div>
</li>
<li class="hidden">
<a class="dashboard-shortcuts-projects" href="/dashboard/projects">Projects
</a></li>
<li class="d-lg-block d-none dropdown">

</li>
<li class="d-none d-lg-block d-xl-block"><a title="Instance Statistics" aria-label="Instance Statistics" data-toggle="tooltip" data-placement="bottom" data-container="body" href="/-/instance_statistics"><svg class="s18"><use xlink:href="/assets/icons-4009ebf96719f129f954d643e65f87152b5e6c4a4917130c4a696beb54af9949.svg#chart"></use></svg>
</a></li>
</ul>

</div>
<div class="navbar-collapse collapse">
<ul class="nav navbar-nav">
<li class="header-new dropdown" data-track-event="click_dropdown" data-track-label="new_dropdown" data-track-value="">
<a class="header-new-dropdown-toggle has-tooltip qa-new-menu-toggle" id="js-onboarding-new-project-link" title="New..." ref="tooltip" aria-label="New..." data-toggle="dropdown" data-placement="bottom" data-container="body" data-display="static" href="/projects/new"><svg class="s16"><use xlink:href="/assets/icons-4009ebf96719f129f954d643e65f87152b5e6c4a4917130c4a696beb54af9949.svg#plus-square"></use></svg>
<svg class="caret-down"><use xlink:href="/assets/icons-4009ebf96719f129f954d643e65f87152b5e6c4a4917130c4a696beb54af9949.svg#angle-down"></use></svg>
</a><div class="dropdown-menu dropdown-menu-right">
<ul>
<li class="dropdown-bold-header">
This project
</li>
<li><a href="/samso008/runhisat/issues/new">New issue</a></li>
<li class="divider"></li>
<li class="dropdown-bold-header">GitLab</li>
<li><a class="qa-global-new-project-link" href="/projects/new">New project</a></li>
<li><a href="/groups/new">New group</a></li>
<li><a class="qa-global-new-snippet-link" href="/snippets/new">New snippet</a></li>
</ul>
</div>
</li>

<li class="nav-item d-none d-sm-none d-md-block m-auto">
<div class="search search-form" data-track-event="activate_form_input" data-track-label="navbar_search" data-track-value="">
<form class="form-inline" action="/search" accept-charset="UTF-8" method="get"><input name="utf8" type="hidden" value="&#x2713;" /><div class="search-input-container">
<div class="search-input-wrap">
<div class="dropdown" data-url="/search/autocomplete">
<input type="search" name="search" id="search" placeholder="Search or jump to…" class="search-input dropdown-menu-toggle no-outline js-search-dashboard-options" spellcheck="false" tabindex="1" autocomplete="off" data-issues-path="/dashboard/issues" data-mr-path="/dashboard/merge_requests" data-qa-selector="search_term_field" aria-label="Search or jump to…" />
<button class="hidden js-dropdown-search-toggle" data-toggle="dropdown" type="button"></button>
<div class="dropdown-menu dropdown-select js-dashboard-search-options">
<div class="dropdown-content"><ul>
<li class="dropdown-menu-empty-item">
<a>
Loading...
</a>
</li>
</ul>
</div><div class="dropdown-loading"><i aria-hidden="true" data-hidden="true" class="fa fa-spinner fa-spin"></i></div>
</div>
<svg class="s16 search-icon"><use xlink:href="/assets/icons-4009ebf96719f129f954d643e65f87152b5e6c4a4917130c4a696beb54af9949.svg#search"></use></svg>
<svg class="s16 clear-icon js-clear-input"><use xlink:href="/assets/icons-4009ebf96719f129f954d643e65f87152b5e6c4a4917130c4a696beb54af9949.svg#close"></use></svg>
</div>
</div>
</div>
<input type="hidden" name="group_id" id="group_id" class="js-search-group-options" />
<input type="hidden" name="project_id" id="search_project_id" value="5460" class="js-search-project-options" data-project-path="runhisat" data-name="runhisat" data-issues-path="/samso008/runhisat/issues" data-mr-path="/samso008/runhisat/merge_requests" data-issues-disabled="false" />
<input type="hidden" name="search_code" id="search_code" value="true" />
<input type="hidden" name="repository_ref" id="repository_ref" value="master" />
<input type="hidden" name="nav_source" id="nav_source" value="navbar" />
<div class="search-autocomplete-opts hide" data-autocomplete-path="/search/autocomplete" data-autocomplete-project-id="5460" data-autocomplete-project-ref="master"></div>
</form></div>

</li>
<li class="nav-item d-inline-block d-sm-none d-md-none">
<a title="Search" aria-label="Search" data-toggle="tooltip" data-placement="bottom" data-container="body" href="/search?project_id=5460"><svg class="s16"><use xlink:href="/assets/icons-4009ebf96719f129f954d643e65f87152b5e6c4a4917130c4a696beb54af9949.svg#search"></use></svg>
</a></li>
<li class="user-counter"><a title="Issues" class="dashboard-shortcuts-issues" aria-label="Issues" data-toggle="tooltip" data-placement="bottom" data-container="body" href="/dashboard/issues?assignee_username=lohr001"><svg class="s16"><use xlink:href="/assets/icons-4009ebf96719f129f954d643e65f87152b5e6c4a4917130c4a696beb54af9949.svg#issues"></use></svg>
<span class="badge badge-pill green-badge hidden issues-count">
0
</span>
</a></li><li class="user-counter"><a title="Merge requests" class="dashboard-shortcuts-merge_requests" aria-label="Merge requests" data-toggle="tooltip" data-placement="bottom" data-container="body" href="/dashboard/merge_requests?assignee_username=lohr001"><svg class="s16"><use xlink:href="/assets/icons-4009ebf96719f129f954d643e65f87152b5e6c4a4917130c4a696beb54af9949.svg#git-merge"></use></svg>
<span class="badge badge-pill hidden merge-requests-count">
0
</span>
</a></li><li class="user-counter"><a title="To-Do List" aria-label="To-Do List" class="shortcuts-todos" data-toggle="tooltip" data-placement="bottom" data-container="body" href="/dashboard/todos"><svg class="s16"><use xlink:href="/assets/icons-4009ebf96719f129f954d643e65f87152b5e6c4a4917130c4a696beb54af9949.svg#todo-done"></use></svg>
<span class="badge badge-pill hidden todos-count">
0
</span>
</a></li><li class="nav-item header-help dropdown">
<a class="header-help-dropdown-toggle" data-toggle="dropdown" href="/help"><svg class="s16"><use xlink:href="/assets/icons-4009ebf96719f129f954d643e65f87152b5e6c4a4917130c4a696beb54af9949.svg#question"></use></svg>
<svg class="caret-down"><use xlink:href="/assets/icons-4009ebf96719f129f954d643e65f87152b5e6c4a4917130c4a696beb54af9949.svg#angle-down"></use></svg>
</a><div class="dropdown-menu dropdown-menu-right">
<ul>
<li>
<a href="/help">Help</a>
</li>
<li>
<a href="https://about.gitlab.com/getting-help/">Support</a>
</li>

<li class="divider"></li>
<li>
<a href="https://about.gitlab.com/submit-feedback">Submit feedback</a>
</li>
<li>
<a target="_blank" class="text-nowrap" href="https://about.gitlab.com/contributing">Contribute to GitLab
</a></li>


</ul>

</div>
</li>
<li class="dropdown header-user nav-item" data-qa-selector="user_menu" data-track-event="click_dropdown" data-track-label="profile_dropdown" data-track-value="">
<a class="header-user-dropdown-toggle" data-toggle="dropdown" href="/lohr001"><img width="23" height="23" class="header-user-avatar qa-user-avatar lazy" data-src="https://secure.gravatar.com/avatar/6e7f9ff27c82118fa1299609f0124e5f?s=46&amp;d=identicon" src="data:image/gif;base64,R0lGODlhAQABAAAAACH5BAEKAAEALAAAAAABAAEAAAICTAEAOw==" />
<svg class="caret-down"><use xlink:href="/assets/icons-4009ebf96719f129f954d643e65f87152b5e6c4a4917130c4a696beb54af9949.svg#angle-down"></use></svg>
</a><div class="dropdown-menu dropdown-menu-right">
<ul>
<li class="current-user">
<div class="user-name bold">
Lohr, Julia
</div>
@lohr001
</li>
<li class="divider"></li>
<li>
<div class="js-set-status-modal-trigger" data-has-status="false"></div>
</li>
<li>
<a class="profile-link" data-user="lohr001" href="/lohr001">Profile</a>
</li>
<li>
<a data-qa-selector="settings_link" href="/profile">Settings</a>
</li>
<li class="divider"></li>
<li>
<a class="sign-out-link" data-qa-selector="sign_out_link" href="/users/sign_out">Sign out</a>
</li>
</ul>

</div>
</li>
</ul>
</div>
<button class="navbar-toggler d-block d-sm-none" type="button">
<span class="sr-only">Toggle navigation</span>
<svg class="s12 more-icon js-navbar-toggle-right"><use xlink:href="/assets/icons-4009ebf96719f129f954d643e65f87152b5e6c4a4917130c4a696beb54af9949.svg#ellipsis_h"></use></svg>
<svg class="s12 close-icon js-navbar-toggle-left"><use xlink:href="/assets/icons-4009ebf96719f129f954d643e65f87152b5e6c4a4917130c4a696beb54af9949.svg#close"></use></svg>
</button>
</div>
</div>
</header>
<div class="js-set-status-modal-wrapper" data-current-emoji="" data-current-message=""></div>

<div class="layout-page page-with-contextual-sidebar">
<div class="nav-sidebar">
<div class="nav-sidebar-inner-scroll">
<div class="context-header">
<a title="runhisat" href="/samso008/runhisat"><div class="avatar-container rect-avatar s40 project-avatar">
<div class="avatar s40 avatar-tile identicon bg1">R</div>
</div>
<div class="sidebar-context-title">
runhisat
</div>
</a></div>
<ul class="sidebar-top-level-items">
<li class="home"><a class="shortcuts-project rspec-project-link" data-qa-selector="project_link" href="/samso008/runhisat"><div class="nav-icon-container">
<svg><use xlink:href="/assets/icons-4009ebf96719f129f954d643e65f87152b5e6c4a4917130c4a696beb54af9949.svg#home"></use></svg>
</div>
<span class="nav-item-name">
Project
</span>
</a><ul class="sidebar-sub-level-items">
<li class="fly-out-top-item"><a href="/samso008/runhisat"><strong class="fly-out-top-item-name">
Project
</strong>
</a></li><li class="divider fly-out-top-item"></li>
<li class=""><a title="Project details" class="shortcuts-project" href="/samso008/runhisat"><span>Details</span>
</a></li><li class=""><a title="Activity" class="shortcuts-project-activity" data-qa-selector="activity_link" href="/samso008/runhisat/activity"><span>Activity</span>
</a></li><li class=""><a title="Releases" class="shortcuts-project-releases" href="/samso008/runhisat/-/releases"><span>Releases</span>
</a></li><li class=""><a title="Cycle Analytics" class="shortcuts-project-cycle-analytics" href="/samso008/runhisat/cycle_analytics"><span>Cycle Analytics</span>
</a></li>
</ul>
</li><li class="active"><a class="shortcuts-tree qa-project-menu-repo" href="/samso008/runhisat/tree/master"><div class="nav-icon-container">
<svg><use xlink:href="/assets/icons-4009ebf96719f129f954d643e65f87152b5e6c4a4917130c4a696beb54af9949.svg#doc-text"></use></svg>
</div>
<span class="nav-item-name" id="js-onboarding-repo-link">
Repository
</span>
</a><ul class="sidebar-sub-level-items">
<li class="fly-out-top-item active"><a href="/samso008/runhisat/tree/master"><strong class="fly-out-top-item-name">
Repository
</strong>
</a></li><li class="divider fly-out-top-item"></li>
<li class="active"><a href="/samso008/runhisat/tree/master">Files
</a></li><li class=""><a id="js-onboarding-commits-link" href="/samso008/runhisat/commits/master">Commits
</a></li><li class=""><a class="qa-branches-link" id="js-onboarding-branches-link" href="/samso008/runhisat/-/branches">Branches
</a></li><li class=""><a href="/samso008/runhisat/-/tags">Tags
</a></li><li class=""><a href="/samso008/runhisat/-/graphs/master">Contributors
</a></li><li class=""><a href="/samso008/runhisat/-/network/master">Graph
</a></li><li class=""><a href="/samso008/runhisat/compare?from=master&amp;to=master">Compare
</a></li><li class=""><a href="/samso008/runhisat/-/graphs/master/charts">Charts
</a></li>
</ul>
</li><li class=""><a class="shortcuts-issues qa-issues-item" href="/samso008/runhisat/issues"><div class="nav-icon-container">
<svg><use xlink:href="/assets/icons-4009ebf96719f129f954d643e65f87152b5e6c4a4917130c4a696beb54af9949.svg#issues"></use></svg>
</div>
<span class="nav-item-name" id="js-onboarding-issues-link">
Issues
</span>
<span class="badge badge-pill count issue_counter">
0
</span>
</a><ul class="sidebar-sub-level-items">
<li class="fly-out-top-item"><a href="/samso008/runhisat/issues"><strong class="fly-out-top-item-name">
Issues
</strong>
<span class="badge badge-pill count issue_counter fly-out-badge">
0
</span>
</a></li><li class="divider fly-out-top-item"></li>
<li class=""><a title="Issues" href="/samso008/runhisat/issues"><span>
List
</span>
</a></li><li class=""><a title="Boards" data-qa-selector="issue_boards_link" href="/samso008/runhisat/-/boards"><span>
Boards
</span>
</a></li><li class=""><a title="Labels" class="qa-labels-link" href="/samso008/runhisat/-/labels"><span>
Labels
</span>
</a></li>
<li class=""><a title="Milestones" class="qa-milestones-link" href="/samso008/runhisat/-/milestones"><span>
Milestones
</span>
</a></li></ul>
</li><li class=""><a class="shortcuts-merge_requests" data-qa-selector="merge_requests_link" href="/samso008/runhisat/merge_requests"><div class="nav-icon-container">
<svg><use xlink:href="/assets/icons-4009ebf96719f129f954d643e65f87152b5e6c4a4917130c4a696beb54af9949.svg#git-merge"></use></svg>
</div>
<span class="nav-item-name" id="js-onboarding-mr-link">
Merge Requests
</span>
<span class="badge badge-pill count merge_counter js-merge-counter">
0
</span>
</a><ul class="sidebar-sub-level-items is-fly-out-only">
<li class="fly-out-top-item"><a href="/samso008/runhisat/merge_requests"><strong class="fly-out-top-item-name">
Merge Requests
</strong>
<span class="badge badge-pill count merge_counter js-merge-counter fly-out-badge">
0
</span>
</a></li></ul>
</li><li class=""><a class="shortcuts-pipelines qa-link-pipelines rspec-link-pipelines" href="/samso008/runhisat/pipelines"><div class="nav-icon-container">
<svg><use xlink:href="/assets/icons-4009ebf96719f129f954d643e65f87152b5e6c4a4917130c4a696beb54af9949.svg#rocket"></use></svg>
</div>
<span class="nav-item-name" id="js-onboarding-pipelines-link">
CI / CD
</span>
</a><ul class="sidebar-sub-level-items">
<li class="fly-out-top-item"><a href="/samso008/runhisat/pipelines"><strong class="fly-out-top-item-name">
CI / CD
</strong>
</a></li><li class="divider fly-out-top-item"></li>
<li class=""><a title="Pipelines" class="shortcuts-pipelines" href="/samso008/runhisat/pipelines"><span>
Pipelines
</span>
</a></li><li class=""><a title="Jobs" class="shortcuts-builds" href="/samso008/runhisat/-/jobs"><span>
Jobs
</span>
</a></li><li class=""><a title="Schedules" class="shortcuts-builds" href="/samso008/runhisat/pipeline_schedules"><span>
Schedules
</span>
</a></li><li class=""><a title="Charts" class="shortcuts-pipelines-charts" href="/samso008/runhisat/pipelines/charts"><span>
Charts
</span>
</a></li></ul>
</li>
<li class=""><a href="/samso008/runhisat/container_registry"><div class="nav-icon-container">
<svg><use xlink:href="/assets/icons-4009ebf96719f129f954d643e65f87152b5e6c4a4917130c4a696beb54af9949.svg#package"></use></svg>
</div>
<span class="nav-item-name">
Packages
</span>
</a><ul class="sidebar-sub-level-items">
<li class="fly-out-top-item"><a href="/samso008/runhisat/container_registry"><strong class="fly-out-top-item-name">
Packages
</strong>
</a></li><li class="divider fly-out-top-item"></li>
<li class=""><a class="shortcuts-container-registry" title="Container Registry" href="/samso008/runhisat/container_registry"><span>Container Registry</span>
</a></li></ul>
</li>
<li class=""><a class="shortcuts-wiki" data-qa-selector="wiki_link" href="/samso008/runhisat/wikis/home"><div class="nav-icon-container">
<svg><use xlink:href="/assets/icons-4009ebf96719f129f954d643e65f87152b5e6c4a4917130c4a696beb54af9949.svg#book"></use></svg>
</div>
<span class="nav-item-name">
Wiki
</span>
</a><ul class="sidebar-sub-level-items is-fly-out-only">
<li class="fly-out-top-item"><a href="/samso008/runhisat/wikis/home"><strong class="fly-out-top-item-name">
Wiki
</strong>
</a></li></ul>
</li><li class=""><a class="shortcuts-snippets" href="/samso008/runhisat/snippets"><div class="nav-icon-container">
<svg><use xlink:href="/assets/icons-4009ebf96719f129f954d643e65f87152b5e6c4a4917130c4a696beb54af9949.svg#snippet"></use></svg>
</div>
<span class="nav-item-name">
Snippets
</span>
</a><ul class="sidebar-sub-level-items is-fly-out-only">
<li class="fly-out-top-item"><a href="/samso008/runhisat/snippets"><strong class="fly-out-top-item-name">
Snippets
</strong>
</a></li></ul>
</li><li class=""><a title="Members" class="shortcuts-tree" href="/samso008/runhisat/-/settings/members"><div class="nav-icon-container">
<svg><use xlink:href="/assets/icons-4009ebf96719f129f954d643e65f87152b5e6c4a4917130c4a696beb54af9949.svg#users"></use></svg>
</div>
<span class="nav-item-name">
Members
</span>
</a><ul class="sidebar-sub-level-items is-fly-out-only">
<li class="fly-out-top-item"><a href="/samso008/runhisat/-/project_members"><strong class="fly-out-top-item-name">
Members
</strong>
</a></li></ul>
</li><a class="toggle-sidebar-button js-toggle-sidebar qa-toggle-sidebar rspec-toggle-sidebar" role="button" title="Toggle sidebar" type="button">
<svg class="icon-angle-double-left"><use xlink:href="/assets/icons-4009ebf96719f129f954d643e65f87152b5e6c4a4917130c4a696beb54af9949.svg#angle-double-left"></use></svg>
<svg class="icon-angle-double-right"><use xlink:href="/assets/icons-4009ebf96719f129f954d643e65f87152b5e6c4a4917130c4a696beb54af9949.svg#angle-double-right"></use></svg>
<span class="collapse-text">Collapse sidebar</span>
</a>
<button name="button" type="button" class="close-nav-button"><svg class="s16"><use xlink:href="/assets/icons-4009ebf96719f129f954d643e65f87152b5e6c4a4917130c4a696beb54af9949.svg#close"></use></svg>
<span class="collapse-text">Close sidebar</span>
</button>
<li class="hidden">
<a title="Activity" class="shortcuts-project-activity" href="/samso008/runhisat/activity"><span>
Activity
</span>
</a></li>
<li class="hidden">
<a title="Network" class="shortcuts-network" href="/samso008/runhisat/-/network/master">Graph
</a></li>
<li class="hidden">
<a title="Charts" class="shortcuts-repository-charts" href="/samso008/runhisat/-/graphs/master/charts">Charts
</a></li>
<li class="hidden">
<a class="shortcuts-new-issue" href="/samso008/runhisat/issues/new">Create a new issue
</a></li>
<li class="hidden">
<a title="Jobs" class="shortcuts-builds" href="/samso008/runhisat/-/jobs">Jobs
</a></li>
<li class="hidden">
<a title="Commits" class="shortcuts-commits" href="/samso008/runhisat/commits/master">Commits
</a></li>
<li class="hidden">
<a title="Issue Boards" class="shortcuts-issue-boards" href="/samso008/runhisat/-/boards">Issue Boards</a>
</li>
</ul>
</div>
</div>

<div class="content-wrapper">

<div class="mobile-overlay"></div>
<div class="alert-wrapper">






<nav class="breadcrumbs container-fluid container-limited" role="navigation">
<div class="breadcrumbs-container">
<button name="button" type="button" class="toggle-mobile-nav"><span class="sr-only">Open sidebar</span>
<i aria-hidden="true" data-hidden="true" class="fa fa-bars"></i>
</button><div class="breadcrumbs-links js-title-container">
<ul class="list-unstyled breadcrumbs-list js-breadcrumbs-list">
<li><a href="/samso008">Samson, Jonathan</a><svg class="s8 breadcrumbs-list-angle"><use xlink:href="/assets/icons-4009ebf96719f129f954d643e65f87152b5e6c4a4917130c4a696beb54af9949.svg#angle-right"></use></svg></li> <li><a href="/samso008/runhisat"><span class="breadcrumb-item-text js-breadcrumb-item-text">runhisat</span></a><svg class="s8 breadcrumbs-list-angle"><use xlink:href="/assets/icons-4009ebf96719f129f954d643e65f87152b5e6c4a4917130c4a696beb54af9949.svg#angle-right"></use></svg></li>

<li>
<h2 class="breadcrumbs-sub-title"><a href="/samso008/runhisat/blob/master/README.md">Repository</a></h2>
</li>
</ul>
</div>

</div>
</nav>

<div class="d-flex"></div>
</div>
<div class="container-fluid container-limited ">
<div class="content" id="content-body">
<div class="flash-container flash-container-page sticky">
</div>

<div class="js-signature-container" data-signatures-path="/samso008/runhisat/commits/d537c6b471ba1277401b19b6c7da405699bac10a/signatures"></div>

<div class="tree-holder" id="tree-holder">
<div class="nav-block">
<div class="tree-ref-container">
<div class="tree-ref-holder">
<form class="project-refs-form" action="/samso008/runhisat/refs/switch" accept-charset="UTF-8" method="get"><input name="utf8" type="hidden" value="&#x2713;" /><input type="hidden" name="destination" id="destination" value="blob" />
<input type="hidden" name="path" id="path" value="README.md" />
<div class="dropdown">
<button class="dropdown-menu-toggle js-project-refs-dropdown qa-branches-select" type="button" data-toggle="dropdown" data-selected="master" data-ref="master" data-refs-url="/samso008/runhisat/refs?sort=updated_desc" data-field-name="ref" data-submit-form-on-click="true" data-visit="true"><span class="dropdown-toggle-text ">master</span><i aria-hidden="true" data-hidden="true" class="fa fa-chevron-down"></i></button>
<div class="dropdown-menu dropdown-menu-paging dropdown-menu-selectable git-revision-dropdown qa-branches-dropdown">
<div class="dropdown-page-one">
<div class="dropdown-title"><span>Switch branch/tag</span><button class="dropdown-title-button dropdown-menu-close" aria-label="Close" type="button"><i aria-hidden="true" data-hidden="true" class="fa fa-times dropdown-menu-close-icon"></i></button></div>
<div class="dropdown-input"><input type="search" id="" class="dropdown-input-field qa-dropdown-input-field" placeholder="Search branches and tags" autocomplete="off" /><i aria-hidden="true" data-hidden="true" class="fa fa-search dropdown-input-search"></i><i aria-hidden="true" data-hidden="true" role="button" class="fa fa-times dropdown-input-clear js-dropdown-input-clear"></i></div>
<div class="dropdown-content"></div>
<div class="dropdown-loading"><i aria-hidden="true" data-hidden="true" class="fa fa-spinner fa-spin"></i></div>
</div>
</div>
</div>
</form>
</div>
<ul class="breadcrumb repo-breadcrumb">
<li class="breadcrumb-item">
<a href="/samso008/runhisat/tree/master">runhisat
</a></li>
<li class="breadcrumb-item">
<a href="/samso008/runhisat/blob/master/README.md"><strong>README.md</strong>
</a></li>
</ul>
</div>
<div class="tree-controls">
<a class="btn shortcuts-find-file" rel="nofollow" href="/samso008/runhisat/find_file/master"><i aria-hidden="true" data-hidden="true" class="fa fa-search"></i>
<span>Find file</span>
</a>
<div class="btn-group" role="group"><a class="btn js-blob-blame-link" href="/samso008/runhisat/blame/master/README.md">Blame</a><a class="btn" href="/samso008/runhisat/commits/master/README.md">History</a><a class="btn js-data-file-blob-permalink-url" href="/samso008/runhisat/blob/d537c6b471ba1277401b19b6c7da405699bac10a/README.md">Permalink</a></div>
</div>
</div>

<div class="info-well d-none d-sm-block">
<div class="well-segment">
<ul class="blob-commit-info">
<li class="commit flex-row js-toggle-container" id="commit-d537c6b4">
<div class="avatar-cell d-none d-sm-block">
<a href="/samso008"><img alt="Samson, Jonathan&#39;s avatar" src="https://secure.gravatar.com/avatar/b056b6f5bf208d731d5e7d3935df7167?s=80&amp;d=identicon" class="avatar s40 d-none d-sm-inline-block" title="Samson, Jonathan" /></a>
</div>
<div class="commit-detail flex-list">
<div class="commit-content qa-commit-content">
<a class="commit-row-message item-title js-onboarding-commit-item" href="/samso008/runhisat/commit/d537c6b471ba1277401b19b6c7da405699bac10a">Changed a reference to runkallisto to runhisat.</a>
<span class="commit-row-message d-inline d-sm-none">
&middot;
d537c6b4
</span>
<div class="committer">
<a class="commit-author-link js-user-link" data-user-id="1252" href="/samso008">Samson, Jonathan</a> authored <time class="js-timeago" title="Aug 29, 2019 5:43am" datetime="2019-08-29T05:43:34Z" data-toggle="tooltip" data-placement="bottom" data-container="body">Aug 29, 2019</time>
</div>

</div>
<div class="commit-actions flex-row">

<div class="js-commit-pipeline-status" data-endpoint="/samso008/runhisat/commit/d537c6b471ba1277401b19b6c7da405699bac10a/pipelines?ref=master"></div>
<div class="commit-sha-group d-none d-sm-flex">
<div class="label label-monospace monospace">
d537c6b4
</div>
<button class="btn btn btn-default" data-toggle="tooltip" data-placement="bottom" data-container="body" data-title="Copy commit SHA to clipboard" data-class="btn btn-default" data-clipboard-text="d537c6b471ba1277401b19b6c7da405699bac10a" type="button" title="Copy commit SHA to clipboard" aria-label="Copy commit SHA to clipboard"><svg><use xlink:href="/assets/icons-4009ebf96719f129f954d643e65f87152b5e6c4a4917130c4a696beb54af9949.svg#duplicate"></use></svg></button>

</div>
</div>
</div>
</li>

</ul>
</div>


</div>
<div class="blob-content-holder" id="blob-content-holder">
<article class="file-holder">
<div class="js-file-title file-title-flex-parent">
<div class="file-header-content">
<i aria-hidden="true" data-hidden="true" class="fa fa-file-text-o fa-fw"></i>
<strong class="file-title-name qa-file-title-name">
README.md
</strong>
<button class="btn btn-clipboard btn-transparent" data-toggle="tooltip" data-placement="bottom" data-container="body" data-class="btn-clipboard btn-transparent" data-title="Copy file path to clipboard" data-clipboard-text="{&quot;text&quot;:&quot;README.md&quot;,&quot;gfm&quot;:&quot;`README.md`&quot;}" type="button" title="Copy file path to clipboard" aria-label="Copy file path to clipboard"><svg><use xlink:href="/assets/icons-4009ebf96719f129f954d643e65f87152b5e6c4a4917130c4a696beb54af9949.svg#duplicate"></use></svg></button>
<small class="mr-1">
7.8 KB
</small>
</div>

<div class="file-actions">
<div class="btn-group js-blob-viewer-switcher" role="group">
<button aria-label="Display source" class="btn btn-default btn-sm js-blob-viewer-switch-btn has-tooltip" data-container="body" data-viewer="simple" title="Display source">
<i aria-hidden="true" data-hidden="true" class="fa fa-code"></i>
</button><button aria-label="Display rendered file" class="btn btn-default btn-sm js-blob-viewer-switch-btn has-tooltip" data-container="body" data-viewer="rich" title="Display rendered file">
<i aria-hidden="true" data-hidden="true" class="fa fa-file-text-o"></i>
</button></div>

<div class="btn-group" role="group"><button class="btn btn btn-sm js-copy-blob-source-btn" data-toggle="tooltip" data-placement="bottom" data-container="body" data-class="btn btn-sm js-copy-blob-source-btn" data-title="Copy source to clipboard" data-clipboard-target=".blob-content[data-blob-id=&#39;0867774b5f7de0fee4aa11d96868dd941272843b&#39;]" type="button" title="Copy source to clipboard" aria-label="Copy source to clipboard"><svg><use xlink:href="/assets/icons-4009ebf96719f129f954d643e65f87152b5e6c4a4917130c4a696beb54af9949.svg#duplicate"></use></svg></button><a class="btn btn-sm has-tooltip" target="_blank" rel="noopener noreferrer" title="Open raw" data-container="body" href="/samso008/runhisat/raw/master/README.md"><i aria-hidden="true" data-hidden="true" class="fa fa-file-code-o"></i></a><a download="README.md" class="btn btn-sm has-tooltip" target="_blank" rel="noopener noreferrer" title="Download" data-container="body" href="/samso008/runhisat/raw/master/README.md?inline=false"><svg><use xlink:href="/assets/icons-4009ebf96719f129f954d643e65f87152b5e6c4a4917130c4a696beb54af9949.svg#download"></use></svg></a></div>
<div class="btn-group" role="group"><button name="button" type="submit" class="btn js-edit-blob  js-edit-blob-link-fork-toggler" data-action="edit" data-fork-path="/samso008/runhisat/-/forks?continue%5Bnotice%5D=You%27re+not+allowed+to+make+changes+to+this+project+directly.+A+fork+of+this+project+has+been+created+that+you+can+make+changes+in%2C+so+you+can+submit+a+merge+request.&amp;continue%5Bnotice_now%5D=You%27re+not+allowed+to+make+changes+to+this+project+directly.+A+fork+of+this+project+is+being+created+that+you+can+make+changes+in%2C+so+you+can+submit+a+merge+request.&amp;continue%5Bto%5D=%2Fsamso008%2Frunhisat%2Fedit%2Fmaster%2FREADME.md&amp;namespace_key=1614">Edit</button><button name="button" type="submit" class="btn btn-default js-edit-blob-link-fork-toggler" data-action="edit" data-fork-path="/samso008/runhisat/-/forks?continue%5Bnotice%5D=You%27re+not+allowed+to+make+changes+to+this+project+directly.+A+fork+of+this+project+has+been+created+that+you+can+make+changes+in%2C+so+you+can+submit+a+merge+request.&amp;continue%5Bnotice_now%5D=You%27re+not+allowed+to+make+changes+to+this+project+directly.+A+fork+of+this+project+is+being+created+that+you+can+make+changes+in%2C+so+you+can+submit+a+merge+request.&amp;continue%5Bto%5D=%2F-%2Fide%2Fproject%2Flohr001%2Frunhisat%2Fedit%2Fmaster%2F-%2FREADME.md&amp;namespace_key=1614">Web IDE</button><button name="button" type="submit" class="btn btn-default js-edit-blob-link-fork-toggler" data-action="replace" data-fork-path="/samso008/runhisat/-/forks?continue%5Bnotice%5D=You%27re+not+allowed+to+make+changes+to+this+project+directly.+A+fork+of+this+project+has+been+created+that+you+can+make+changes+in%2C+so+you+can+submit+a+merge+request.+Try+to+replace+this+file+again.&amp;continue%5Bnotice_now%5D=You%27re+not+allowed+to+make+changes+to+this+project+directly.+A+fork+of+this+project+is+being+created+that+you+can+make+changes+in%2C+so+you+can+submit+a+merge+request.&amp;continue%5Bto%5D=%2Fsamso008%2Frunhisat%2Fblob%2Fmaster%2FREADME.md&amp;namespace_key=1614">Replace</button><button name="button" type="submit" class="btn btn-remove js-edit-blob-link-fork-toggler" data-action="delete" data-fork-path="/samso008/runhisat/-/forks?continue%5Bnotice%5D=You%27re+not+allowed+to+make+changes+to+this+project+directly.+A+fork+of+this+project+has+been+created+that+you+can+make+changes+in%2C+so+you+can+submit+a+merge+request.+Try+to+delete+this+file+again.&amp;continue%5Bnotice_now%5D=You%27re+not+allowed+to+make+changes+to+this+project+directly.+A+fork+of+this+project+is+being+created+that+you+can+make+changes+in%2C+so+you+can+submit+a+merge+request.&amp;continue%5Bto%5D=%2Fsamso008%2Frunhisat%2Fblob%2Fmaster%2FREADME.md&amp;namespace_key=1614">Delete</button></div>
</div>
</div>
<div class="js-file-fork-suggestion-section file-fork-suggestion hidden">
<span class="file-fork-suggestion-note">
You're not allowed to
<span class="js-file-fork-suggestion-section-action">
edit
</span>
files in this project directly. Please fork this project,
make your changes there, and submit a merge request.
</span>
<a class="js-fork-suggestion-button btn btn-grouped btn-inverted btn-success" rel="nofollow" data-method="post" href="/samso008/runhisat/blob/master/README.md">Fork</a>
<button class="js-cancel-fork-suggestion-button btn btn-grouped" type="button">
Cancel
</button>
</div>



<div class="blob-viewer hidden" data-type="simple" data-url="/samso008/runhisat/blob/master/README.md?format=json&amp;viewer=simple">
<div class="text-center prepend-top-default append-bottom-default">
<i aria-hidden="true" aria-label="Loading content…" class="fa fa-spinner fa-spin fa-2x qa-spinner"></i>
</div>

</div>

<div class="blob-viewer" data-rich-type="markup" data-type="rich" data-url="/samso008/runhisat/blob/master/README.md?format=json&amp;viewer=rich">
<div class="text-center prepend-top-default append-bottom-default">
<i aria-hidden="true" aria-label="Loading content…" class="fa fa-spinner fa-spin fa-2x qa-spinner"></i>
</div>

</div>


</article>
</div>

<div class="modal" id="modal-upload-blob">
<div class="modal-dialog modal-lg">
<div class="modal-content">
<div class="modal-header">
<h3 class="page-title">Replace README.md</h3>
<button aria-label="Close" class="close" data-dismiss="modal" type="button">
<span aria-hidden="true">&times;</span>
</button>
</div>
<div class="modal-body">
<form class="js-quick-submit js-upload-blob-form" data-method="put" action="/samso008/runhisat/update/master/README.md" accept-charset="UTF-8" method="post"><input name="utf8" type="hidden" value="&#x2713;" /><input type="hidden" name="_method" value="put" /><input type="hidden" name="authenticity_token" value="HnM2Uq3+rTH7fALS8D3JPN9rWMWTGxXnRRLrjb33Ac5U67SbBhAksWYzs2Pk5f2Q3sIXB4hrggJSG4sCyNlD/A==" /><div class="dropzone">
<div class="dropzone-previews blob-upload-dropzone-previews">
<p class="dz-message light">
Attach a file by drag &amp; drop or <a class="markdown-selector" href="#">click to upload</a>
</p>
</div>
</div>
<br>
<div class="dropzone-alerts alert alert-danger data" style="display:none"></div>
<div class="form-group row commit_message-group">
<label class="col-form-label col-sm-2" for="commit_message-dcc727719582a1a45f4b8207bdd04f8e">Commit message
</label><div class="col-sm-10">
<div class="commit-message-container">
<div class="max-width-marker"></div>
<textarea name="commit_message" id="commit_message-dcc727719582a1a45f4b8207bdd04f8e" class="form-control js-commit-message" placeholder="Replace README.md" required="required" rows="3">
Replace README.md</textarea>
</div>
</div>
</div>

<input type="hidden" name="branch_name" id="branch_name" />
<input type="hidden" name="create_merge_request" id="create_merge_request" value="1" />
<input type="hidden" name="original_branch" id="original_branch" value="master" class="js-original-branch" />

<div class="form-actions">
<button name="button" type="button" class="btn btn-success btn-upload-file" id="submit-all"><i aria-hidden="true" data-hidden="true" class="fa fa-spin fa-spinner js-loading-icon hidden"></i>
Replace file
</button><a class="btn btn-cancel" data-dismiss="modal" href="#">Cancel</a>
<div class="inline prepend-left-10">
A new branch will be created in your fork and a new merge request will be started.
</div>

</div>
</form></div>
</div>
</div>
</div>

</div>

</div>
</div>
</div>
</div>




</body>
</html>

